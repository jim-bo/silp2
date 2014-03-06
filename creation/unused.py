
def _extract_sam(output):
    ''' extracts output form SAM'''

    # extract unique to file, save multimap annotations.
    for line in iter(output,''):

        # skip header.
        if line[0] == '@': continue

        # split.
        tokens = line.strip().split()

        # check for no align.
        if tokens[2] == '*':
            continue 

        # skip low scores.
        if int(tokens[4]) < 10:
            yield False, line
            continue

        # look for XS tag.
        if line.count("XS:i:") == 0:
            
            # no chance of duplicate.
            yield True, line
            
        else:
            
            # find XS and AS.
            xs = None
            zs = None
            for x in tokens[11:-1]:
                if x.count("XS") > 0:
                    xs = int(x.replace("XS:i:",""))
                if x.count("AS") > 0:
                    zs = int(x.replace("AS:i:",""))

            # skip if non-comformist
            if xs == None or zs == None:
                continue
                
            # yield good or not.
            if xs == zs:
                yield False, line
            else:
                yield True, line


def _sam_gen(map1, map2, key_size):
    '''yields the SAM name and the line index'''

    # loop till end of file.
    line1 = map1.readline()
    line2 = map2.readline()
    pos1 = 0
    pos2 = 0
    while line1 != '' and line2 != '':

        # process it.
        tok1 = line1.strip().split()
        tok2 = line2.strip().split()

        # remove to key.
        if key_size != 0:
            key1 = tok1[0][0:-key_size]
            key2 = tok2[0][0:-key_size]
        else:
            key1 = tok1[0]
            key2 = tok2[0]

        # yield the name and line number.
        yield (pos1, key1), (pos2, key2)

        # update info.
        pos1 += len(line1)
        pos2 += len(line2)
        line1 = map1.readline()
        line2 = map2.readline()

def _pair_gen(hitlist1, hitlist2):
    ''' does an in-order walk to find pairs '''

    # loop till each list is empty.
    while len(hitlist1) > 0 and len(hitlist2) > 0:

        # peek for a match.
        if hitlist1[-1][1] == hitlist2[-1][1]:

            # yield it.
            yield hitlist1[-1], hitlist2[-1]

            # change left.
            hitlist1.pop()

        else:

            # pop smaller.
            if hitlist1[-1][1] < hitlist2[-1][1]:
                hitlist1.pop()
            else:
                hitlist2.pop()

def _numpy_unique(srt1):
    """ return unique subset"""

    # create mask.
    good = np.zeros(srt1.shape[0], dtype=np.bool)
    good[:] = False

    # iterate over non-boundry cases.
    for i in range(1, srt1.shape[0]-1):

        # must not match its neighbors.
        if srt1[i-1] != srt1[i] and srt1[i+1] != srt1[i]:
            good[i] = True

    # check the first one.
    if srt1[0] != srt1[1]:
        good[0] = True

    # check the last one.
    if srt1[-1] != srt1[-2]:
        good[-1] = True

    # return the subset slice.
    return srt1[good]
    
def _write_valid(sam_in_1, id1, valid, sam_out_1):
    """ writes entries from valid set"""

    # open output.
    fout1 = open(sam_out_1, "wb")
    fin1 = open(sam_in_1, "rb")

    # generator of pairs.
    idx = 0
    for line in fin1:
        
        # operate.
        if id1[idx]['name'] in valid:
            fout1.write(line)

        # udpate
        idx += 1

    # close em.
    fout1.close()
    fin1.close()

def _uq_mmap(mmap_name, name_size, mode, length=None):
    """ returns pointer to mapped file """
    
    if length == None:
        return np.memmap(mmap_name, dtype='S%d' % name_size, mode=mode)
    else:
        return np.memmap(mmap_name, dtype='S%d' % name_size, mode=mode, shape=(length,))
        


def pair_sam2(sam_in_1, sam_in_2, sam_out_1, sam_out_2, key_size, base_dir):
    """ pairs SAM files """

    # memory map files.
    mmap_id1 = '%s/id1.dat' % base_dir
    mmap_id2 = '%s/id2.dat' % base_dir
    mmap_uq1 = '%s/uq1.dat' % base_dir
    mmap_uq2 = '%s/uq2.dat' % base_dir

    # build name arrays.
    if os.path.isfile(mmap_id1) == False:
        logging.info("extracting name array 1")
        line_cnt1 = _extract_names(sam_in_1, mmap_id1, key_size)
    
    if os.path.isfile(mmap_id2) == False:
        logging.info("extracting name array 2")
        line_cnt2 = _extract_names(sam_in_2, mmap_id2,  key_size)
        
    # compute name sizes.
    name_size1 = _name_size(sam_in_1)
    name_size2 = _name_size(sam_in_2)
        
    # compute uniques for first.
    if os.path.isfile(mmap_uq1) == False:
        logging.info("loading name array 1")
        tmp = _name_mmap(mmap_id1, name_size1, 'c', length=None)
        srt1 = tmp.copy()
        del tmp
        
        logging.info("sorting name array 1")
        srt1.sort(order=['name'])
        
        logging.info("find unique 1")
        uq1 = _numpy_unique(srt1)
        del srt1
        
        logging.info("serializing unique 1")
        tmp = _uq_mmap(mmap_uq1, name_size1, 'write', length=uq1.shape[0])
        tmp[:] = uq1[:]['name']
        del tmp
        del uq1
    
    # compute uniques for second.
    if os.path.isfile(mmap_uq2) == False:
        logging.info("loading name array 2")
        tmp = _name_mmap(mmap_id2, name_size2, 'c', length=None)
        srt2 = tmp.copy()
        
        logging.info("sorting name array 2")
        srt2.sort(order=['name'])
        
        logging.info("find unique 2")
        uq2 = _numpy_unique(srt2)
        del srt2
        
        logging.info("serializing unique 2")
        tmp = _uq_mmap(mmap_uq2, name_size2, 'write', length=uq2.shape[0])
        tmp[:] = uq2[:]['name']
        del tmp
        del uq2

    # compute the intersection of unique.
    logging.info("finding intersection")
    tmp = _uq_mmap(mmap_uq1, name_size1, 'c')
    uq1 = tmp.copy()
    del tmp
    tmp = _uq_mmap(mmap_uq2, name_size2, 'c')
    uq2 = tmp.copy()
    del tmp
    
    valid_list = np.intersect1d(uq1, uq2, assume_unique=True)
    del uq1
    del uq2
    
    # sanity check.
    assert len(valid_list) != 0, 'cant have no valid stuff'

    # create a set.
    logging.info("cast to set")
    valid = set(list(valid_list))

    # write the entries.
    logging.info("writing aligned SAM file 1")
    tmp = _name_mmap(mmap_id1, name_size1, 'r', length=None)
    id1 = tmp.copy()
    _write_valid(sam_in_1, id1, valid, sam_out_1)
    del id1
    del tmp
    
    logging.info("writing aligned SAM file 2")
    tmp = _name_mmap(mmap_id2, name_size2, 'r', length=None)
    id2 = tmp.copy()
    _write_valid(sam_in_2, id2, valid, sam_out_2)
    del id2
    del tmp
    
    logging.info("done")


def pair_sam(sam_in_1, sam_in_2, sam_out_1, sam_out_2, key_size):
    """ pairs SAM files """

    # memory map the SAM files1.
    fin1 = open(sam_in_1, "r+")
    fin2 = open(sam_in_2, "r+")

    map1 = mmap.mmap(fin1.fileno(), 0, access=mmap.ACCESS_COPY)
    map2 = mmap.mmap(fin2.fileno(), 0, access=mmap.ACCESS_COPY)

    # create lists from data.
    hitlist1 = list()
    hitlist2 = list()
    for p1, p2 in _sam_gen(map1, map2, key_size):
        hitlist1.append(p1)
        hitlist2.append(p2)

    # seek files bake to begining.
    map1.seek(0)
    map2.seek(0)

    # sort lists by name, reverse so we can pop from end.
    hitlist1.sort(key=itemgetter(1), reverse=True)
    hitlist2.sort(key=itemgetter(1), reverse=True)

    # open output files.
    fout1 = open(sam_out_1, "wb")
    fout2 = open(sam_out_2, "wb")

    # generator of pairs.
    for p1, p2 in _pair_gen(hitlist1, hitlist2):

        # load sam info from map.
        map1.seek(p1[0])
        map2.seek(p2[0])

        # write out info.
        fout1.write(map1.readline())
        fout2.write(map2.readline())

    # close output files.
    fout1.close()
    fout2.close()

    # close memmory mapped files.
    map1.close()
    map2.close()

    fin1.close()
    fin2.close()
