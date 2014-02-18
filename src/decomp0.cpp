/*
 * decompose.cpp
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

// header
#include "decompose.h"

// namespaces
using namespace std;
using namespace ogdf;

// I/O buffer.
char BUFFER [512];

/*
 * helper functions
 */

// breaks string at tab.
void tokenize(string txt, vector<string> & result){
	
	// clear result.
	result.clear();
	
	// define seperator.
	boost::char_separator<char> sep("\t");
	
	// tokenize the line.
	boost::tokenizer<boost::char_separator<char> > tokens(txt, sep);
	
	// add tokens to vector.
	BOOST_FOREACH(string t, tokens){
		result.push_back(t);
	}
}

Graph load_graph(const char * file_path){
	Graph G;
	string line;
	vector<string> tokens;
	ifstream fin(file_path);
	if( fin.is_open() == true ){
		
		// grab the header.
		getline(fin, line);
		tokenize(line, tokens);
		
		// parse.
		int num_nodes = atoi(tokens[0].c_str());
		int num_edges = atoi(tokens[1].c_str());
		
		// add the nodes.
		Array<node> nlist(0, num_nodes);		
		for(int i=0; i<num_nodes; i++){
			nlist[i] = G.newNode();
		}
		
		// add the edges.
		int p, q;
		for(int i=0; i<num_edges; i++){
			getline(fin, line);
			tokenize(line, tokens);	
			p = atoi(tokens[0].c_str());
			q = atoi(tokens[1].c_str());
			G.newEdge(nlist[p], nlist[q]);
		}
	} else {
		WRITE_ERR("bad input file\n");
		exit(1);
	}
	
	// close the file.
	fin.close();
	
	// return the graph.
	return G;
}

int main(int argc, char* argv[]) {

	// arguments.
	const char * in_file = argv[1];
	const char * out_file = argv[2];
	
	// read the graph.
	WRITE_OUT("loading graph\n");
	Graph G = load_graph(in_file);
	
	// sanity check.
	if( isConnected(G) == false ){
		WRITE_ERR("not connected?\n");
		exit(1);
	}
	
	// perform decomposition.
	WRITE_OUT("performing decomposition\n");
	BCTree D(G, false);
	
	sprintf(BUFFER, "Bnodes: %d\nCnode: %d\n", D.numberOfBComps(), D.numberOfCComps());
	WRITE_OUT(BUFFER);
	
	// get the tree.
	const Graph &T = D.bcTree();
	
	// open output file.
	FILE * fout = fopen(out_file, "w");
	
	// identify root of bcTree
	node root = NULL;
	node n;
	edge e;
	int cnt;
	forall_nodes(n, T){
		
		// count outgoing edges.
		cnt = 0;
		forall_adj_edges(e, n) {
			if(e->source() != n) continue;
			cnt += 1;
		}
		
		// die once we find it.		
		if( cnt == 0 ){
			root = n;
			break;
		}
	}
	
	// sanity.
	if( root == NULL){
		WRITE_ERR("couldn't find a root in bcTREE\n");
		exit(1);
	}

	// write the node entries and sets.
	WRITE_OUT("writing components\n");
	map<int, int> compmap;
	set<int> comp;
	set<int>::iterator cit;
	SListIterator<edge> eit;
	int cidx = 0;
	forall_nodes(n, T){
		if( D.typeOfBNode(n) == D.BComp ){
			
			// get a list of edges.
			SList<edge> cedges = D.hEdges(n);
			
			// build a set of original nodes.
			comp.clear();
			for(eit=cedges.begin(); eit.valid(); ++eit){
				e = *eit;
				comp.insert(D.original(e->source())->index());
				comp.insert(D.original(e->target())->index());
			}
			
			// assign id.
			compmap[n->index()] = cidx;
			cidx++;
			
			// write out id.
			fprintf(fout, "N\t%d", compmap[n->index()]);
			
			// write out sets.
			for(cit=comp.begin(); cit!=comp.end(); cit++){
				fprintf(fout, "\t%d", *cit);
			}
			fprintf(fout, "\n");
		}		
	}
	
	// free up memory.
	comp.clear();
	
	// write the cut edges.
	WRITE_OUT("writing cuts\n");
	set<int> parents, children;
	set<int>::iterator pit;
	edge e1, e2;
	node s, t, c;
	forall_nodes(n, T){
		if( D.typeOfBNode(n) == D.CComp ){
					
			// build set of parents.
			parents.clear();
			forall_adj_edges(e1, n) {
				if(e1->source() != n) continue;
				parents.insert(compmap[e1->target()->index()]);
			}
			
			// build set of children.
			children.clear();
			forall_adj_edges(e1, n) {
				if(e1->target() != n) continue;
				children.insert(compmap[e1->source()->index()]);
			}
					
			// identify the cut-node.
			forall_adj_edges(e1, n) {
				if(e1->source() != n) continue;
				c = D.original(D.cutVertex(n, e1->source()));
				break;
			}	
			
			// write all combinations.
			for(pit=parents.begin(); pit!=parents.end(); pit++){
				for(cit=children.begin(); cit!=children.end(); cit++){
					fprintf(fout, "E\t%d\t%d\t%d\n", *pit, *cit, c->index());
				}
			}
		}
	}	
	
	// close output.
	fclose(fout);
}
