#include "Graph.h"
#include <algorithm>



/* constructor, taking a simple graph class */
MVGraph::MVGraph(graph& g) {
	matchnum = 0;
	bridgenum = 0;
	todonum = 0;
	
	/* simply move edges from the other graph */
	edges = g.edges;
	g.edges = 0;
	if(!g.edges_owned) v_edges.swap(g.edges_vect);
	edges_size = g.nedges;
	
	/* allocate nodes, copy important parameters */
	nodes.resize(g.nnodes);
	for(nodeid i=0;i<g.nnodes;i++) {
		nodes[i].edges_idx = g.idx[i];
		nodes[i].degree = g.outdeg[i];
	}
	/* free all memory associated with the previous graph class */
	g.clear();
	
	levels.reserve(nodes.size()/2+1);
	bridges.reserve(nodes.size()/2+1);	
}



void MVGraph::add_to_level(unsigned int level, nodeid node) {
	if(level >= levels.size()) levels.resize(level+1);
	levels[level].push_back(node);
  //~ add_to_list(get(g->levels,level),node);
	todonum++;
}


void MVGraph::add_to_bridges(unsigned int level, nodeid n1, nodeid n2) {
	if(level >= bridges.size()) bridges.resize(level+1);
	bridges[level].push_back(MVBridge(n1,n2));
	bridgenum++;;
}


void MVGraph::reset() {
	for(auto& v : levels) v.clear();
	for(auto& v : bridges) v.clear();
	bridgenum = 0;
	todonum = 0;
	for(nodeid i=0;i<nodes.size();i++) {
		nodes[i].reset();
		if(nodes[i].match == UNMATCHED) {
			add_to_level(0,i);
			nodes[i].set_min_level(0);
		}
	}
}

void MVGraph::greedy_init() {
	for(nodeid j = 0; j < nodes.size(); j++) {
		MVNode& n = nodes[j];
		if(n.match == UNMATCHED) {
			for(unsigned int k=0;k<n.degree;k++) {
				nodeid i = edges[n.edges_idx + k];
				if(nodes[i].match == UNMATCHED) {
					n.match = i;
					nodes[i].match = j; // note: j is size_t
					matchnum++;
					break;
				}
			}
		}
	}
}



void MVGraph::write_matches(FILE* f) const {
	/* write out matches, but only one way for each edge */
	for(nodeid i=0;i<nodes.size();i++) 
		if(nodes[i].match != UNMATCHED && nodes[i].match > i)
			fprintf(f,"%u\t%u\n",i,nodes[i].match);
}

void MVGraph::write_matches(FILE* f, std::function<unsigned int(nodeid)> n) const {
	/* write out matches, but only one way for each edge */
	for(nodeid i=0;i<nodes.size();i++) 
		if(nodes[i].match != UNMATCHED && nodes[i].match > i)
			fprintf(f,"%u\t%u\n",n(i),n(nodes[i].match));
}
	

