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
	if( isBiconnected(G) == false ){
		WRITE_ERR("not biconnected?\n");
		exit(1);
	}
	
	// perform decomposition.
	WRITE_OUT("performing decomposition\n");
	StaticSPQRTree D(G);
	
	sprintf(BUFFER, "Snodes: %d\nPnode: %d\nRnode: %d\n", D.numberOfSNodes(), D.numberOfPNodes(), D.numberOfRNodes());
	WRITE_OUT(BUFFER);
	
	// get the tree.
	const Graph &T = D.tree();
	
	// open output file.
	FILE * fout = fopen(out_file, "w");
	
	// write the node entries and sets.
	WRITE_OUT("writing components\n");
	node n, r;
	map<int, int> compmap;
	set<int> comp;
	set<int>::iterator cit;
	int cidx = 0;
	
	forall_nodes(n, T){
	
		// get skeleton and its graph.
		const Skeleton &sk = D.skeleton(n);
		const Graph &g = sk.getGraph();
	
		// add nodes from skeleton graph.
		comp.clear();
		forall_nodes(r, g){
			comp.insert(sk.original(r)->index());
		}
		
		// map component to idx.
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
	
	// write out the cuts.
	edge e, se;
	node s, t;
	int ids, idt, idc1, idc2;
	forall_edges(e, T){
		
		// simplify.
		s = e->source();
		t = e->target();
		ids = compmap[s->index()];
		idt = compmap[t->index()];
		
		// get source skeleton and its graph.
		const Skeleton &sk = D.skeleton(s);
		const Graph &g = sk.getGraph();	
		
		// loop over skeleton edges.
		forall_edges(se, g){
			
			// look for virtual whos twin is t.
			if(sk.twinTreeNode(se) == t ){
				
				// simplify.
				idc1 = sk.original(se->source())->index();
				idc2 = sk.original(se->target())->index();
				
				// write out info.
				fprintf(fout, "E\t%d\t%d\t%d\t%d\n", ids, idt, idc1, idc2);
			}
		}
	}
	
	// close output.
	fclose(fout);
}
