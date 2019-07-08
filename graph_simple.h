/*  -*- C++ -*-
 * graph_simple.h -- simple graph implementation storing it as a list of edges
 * 	helper to load graph used for symmetric maximum matching
 * 
 * Copyright 2018 Daniel Kondor <kondor.dani@gmail.com>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * Neither the name of the  nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * 
 */
 
#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include "tsv_iterator.h"
#include "iterator_zip.h"


/* 
 * main class storing a graph and implementing search for connected components
 * and maximum matching */
class graph {
	public:
		graph():nnodes(0),nedges(0),edges(0),outdeg(0),idx(0),edges_size(0),nodes_size(0),edges_owned(true) {  }
		void clear() {
			if(edges_owned) if(edges) free(edges);
			edges_vect.clear(); edges = 0;
			if(outdeg) free(outdeg); if(idx) free(idx); outdeg = 0; idx = 0;
			nnodes = 0; nedges = 0; edges_size = 0; nodes_size = 0;
		}
		~graph() { clear(); }
	
		
		/* There are four ways to create the graph, main difference is the format of input and thus
		 * the memory requirement.
		 * Storing the graph requires 4 bytes / edge + 12 bytes / node. The following functions are
		 * provided for creating the graph:
		 * 
		 * 	1. From a partitioned (sorted) input, supplied as forward iterators of std::pair<unsigned int,unsigned int>.
		 * 		In this case, all edges of a node in the input should be grouped together at least
		 * 		(or already sorted).
		 * 		Ideally, the iterators supplied need not belong to a container, but should be
		 * 		generated / read from a file on the fly. In this case, there is no additional
		 * 		memory requirement over storing the graph. */
		template<class it, class sent>
		int create_graph_partitioned(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map);
		 
		 /* 2. From general input supplied by random access iterators (to std::pair<unsigned int,unsigned int>).
		 * 		In this case, these iterators are used to sort the input. In this case, the additional
		 * 		memory requirement is that of the container of the iterators given, typically
		 * 		8 bytes / edge. */
		template<class it, class sent>
		int create_graph_in_place_sort(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map);
		 
		 /*  3. From general input supplied by a pair of forward iterators;
		 * 		in this case, all edges are copied into temporary arrays and sorted there, so again
		 * 		the data should be generated / read by the iterators on the fly;
		 * 		in this case, there is an additional 4 bytes / edge memory requirement during creating
		 * 		the graph */
		template<class it, class sent>
		int create_graph_copy_sort(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map);
		  
		 /*	4. from general input supplied as std::vectors / C arrays:
		 * 		edges should be supplied as two separate arrays or vector (i.e. edge i points from
		 * 		e1[i] -> e2[i]); this instance will take over the 2nd array (endpoints of the edges)
		 * 		in this case, there is no additional memory requirement (but creating those arrays
		 * 		will already have had 8 bytes / edge memory requirement, of which 4 bytes / edge can
		 * 		be freed after having created the graph */
		int create_graph_arrays(unsigned int* e1, unsigned int* e2, size_t e_size, std::unordered_map<unsigned int, unsigned int>* ids_map);
		int create_graph_vectors(std::vector<unsigned int>& e1, std::vector<unsigned int>& e2, std::unordered_map<unsigned int, unsigned int>* ids_map);
		
		
		/* Helper function for the above to handle reading the graph from a file.
		 * If partitioned == true, it expects already partitioned input; otherwise
		 * it sorts the input.
		 * If ids != null, it replaces ids in the edges to form a continuous
		 * range 0...N-1, and the real ids from the file are stored in the supplied
		 * vector */
		int read_graph(FILE* f, bool partitioned, std::vector<unsigned int>* ids);
		
		/* make the graph symmetric, i.e. for all a->b edges, make sure the reverse, b->a edge is present as well
		 * returns 0 on success, 1 on error */
		int make_symmetric();
		
		/* write graph to the given file; the user can specify a function to map
		 * the continuous ids (0...N-1) to any real node IDs (e.g. the mapping
		 * created when reading the graph) */
	/*	void write_graph(FILE* f, std::function<unsigned int(unsigned int)> ids = [](const unsigned int& x) {return x;}) {
			if(!f) return;
			for(auto it = edges_begin(); it != edges_end(); ++it) fprintf(f,"%u\t%u\n",ids(it->first),ids(it->second));
		}*/
		/* same, but use a std::vector as the ID mapping 
		void write_graph(FILE* f, const std::vector<unsigned int>& ids) {
			write_graph(f,[ids](const unsigned int& x){return ids[x];});
		}*/
		
		friend class MVGraph;
		
		unsigned int get_nodes_size() const { return nnodes; }
		uint64_t get_edges_size() const { return nedges; }
	
	protected:
		unsigned int nnodes;
		uint64_t nedges;
		unsigned int* edges;
		unsigned int* outdeg;
		uint64_t* idx; // node i will have edges in [ edges[idx[i]] ; edges[idx[i]+outdeg[i]] ) -- note: outdeg could be omitted and only idx used, but this is more convenient
	
		const size_t EDGES_GROW = 33554432; // allocate memory for edges in 128MB chunks
		const size_t NODES_GROW = 262144; // allocate memory for nodes in 1M and 2M chunks
		
		/* size of the edges and node arrays, and functions to grow them */	
		size_t edges_size;
		size_t nodes_size;
		
		std::vector<unsigned int> edges_vect; /* optionally, the edges are stored in an std::vector<int> received from
			the caller; in this case edges == edges_vect.data() and should not be freed later */
		bool edges_owned; /* this flag determines whether the edges array should be freed by the destructor */
		
		int grow_nodes() {
			size_t new_size = nodes_size + NODES_GROW;
			unsigned int* tmp = (unsigned int*)realloc(outdeg,new_size*sizeof(unsigned int));
			if(!tmp) return 1;
			outdeg = tmp;
			uint64_t* tmp2 = (uint64_t*)realloc(idx,new_size*sizeof(uint64_t));
			if(!tmp2) return 1;
			idx = tmp2;
			nodes_size = new_size;
			return 0;
		}
		
		static int grow_edges_s(unsigned int** e, size_t* s, size_t grow) {
			size_t new_size = *s + grow;
			unsigned int* tmp = (unsigned int*)realloc(*e,new_size*sizeof(unsigned int));
			if(!tmp) return 1;
			*e = tmp;
			*s = new_size;
			return 0;
		}
		
		int grow_edges() {
			return grow_edges_s(&edges,&edges_size,EDGES_GROW);
		}
		
		/* common helper interface used internally doing all the work expecting sorted / partitioned input
		 * it and sent should be zip iterators pointing to the data or iterators of std::pair with the input.
		 * If ids_map is not null, ids are replaced and the map is filled with the mapping;
		 * otherwise, ids are allocated up to the maximum of the id found in the input */
		template<class it, class sent>
		int create_graph_sorted(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map, bool copy_out_edges);
	
		/* 3-4. common interface for creating the edges after copying the out-edges to the edges array and using a separate in-edges array */
		int create_graph(unsigned int* in_edges, std::unordered_map<unsigned int,unsigned int>* ids_map);
		
		unsigned int real_deg(unsigned int n) const;
};



/* templated functions for creating the graph */

template<class it, class sent>
int graph::create_graph_sorted(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map, bool copy_out_edges) {
	if(grow_nodes()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
	if(copy_out_edges) if(grow_edges()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
	if(ids_map) ids_map->clear();
	
	/* 2. count degrees, replace ids if needed */
	unsigned int last_id = 0;
	unsigned int deg = 0;
	unsigned int i = 0; // node index
	uint64_t j = 0; // edge index
	bool first = true;
	unsigned int max_out_edge = 0;
	bool out_edges_sorted = true;
	
	for(;e != end;++e) {
		/* note: potentially using zip_iterators, so cannot use e->first;
		 * but can use any other iterator returning an std::pair<unsigned int, unsigned int> */
		unsigned int id1 = (*e).first;
		unsigned int id2 = (*e).second;
		
		if(id1 != last_id || first) {
			if(!first) {
				outdeg[i] = deg;
				i++;
				if(i == nodes_size) if(grow_nodes()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
			}
			first = false;
			
			idx[i] = j;
			if(ids_map) {
				auto it1 = ids_map->find(id1);
				if(it1 != ids_map->end()) {
					fprintf(stderr,"graph::create_graph_sorted(): edge origin %u appears in more than one place (input not paritioned)!\n",id1);
					return 1;
				}
				ids_map->insert(std::make_pair(id1,i));
			}
			else {
				if(id1 < last_id) {
					fprintf(stderr,"graph::create_graph_sorted(): input is not sorted at edge %u -- %u!\n",id1,id2);
					return 1;
				}
				while(i < id1) {
					outdeg[i] = 0;
					i++;
					if(i == nodes_size) if(grow_nodes()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
					idx[i] = j;
				}
			}
			deg = 0;
			last_id = id1;
		}
		else {
			if(id2 < edges[j-1]) out_edges_sorted = false;
			else if(id2 == edges[j-1] && copy_out_edges) continue; /* note: omit duplicate edges here as well */
		}
		if(!ids_map) if(id2 > max_out_edge) max_out_edge = id2;
		deg++;
		if(copy_out_edges) {
			if(j == edges_size) if(grow_edges()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
			edges[j] = id2;
		}
		j++;
	}
	outdeg[i] = deg;
	
	nedges = j;
	
	/* replace ids in out edges if using ids_map */
	if(ids_map) {
		for(j=0;j<nedges;j++) {
			auto it2 = ids_map->find(edges[j]);
			if(it2 == ids_map->end()) {
				i++;
				if(i == nodes_size) if(grow_nodes()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
				ids_map->insert(std::make_pair(edges[j],i));
				idx[i] = nedges;
				outdeg[i] = 0;
				edges[j] = i;
			}
			else edges[j] = it2->second;
		}
		i++;
	}
	else for(i++;i<max_out_edge;i++) {
		if(i >= nodes_size) if(grow_nodes()) { fprintf(stderr,"graph::create_graph_sorted(): could not allocate memory!\n"); return 1; }
		idx[i] = nedges;
		outdeg[i] = 0;
	}
	
	/* ensure that out edges are still sorted (required by make_symmetric())
	 * required if we found that the out edges were not sorted for any node
	 * or if we replaced ids in out edges, as in that case, we cannot be sure
	 * about the order there */
	if(!out_edges_sorted || ids_map) for(j=0;j<i;j++) if(outdeg[j] > 1) {
		unsigned int* start = edges + idx[j];
		unsigned int* end = edges + idx[j] + outdeg[j];
		std::sort(start,end);
	}
	
	nnodes = i;
	return 0;
}
		

template<class it, class sent>
int graph::create_graph_partitioned(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map) {
	clear();
	return create_graph_sorted(e,end,ids_map,true);
}
		

/* 2. Create graph with sorting the edges using the supplied iterators (which should be random access
 * 	iterators to std::pair<unsigned int, unsigned int>)
 * note: this function modifies the elements referred to by the iterators (by sorting and possible relabeling IDs) */
template<class it, class sent>
int graph::create_graph_in_place_sort(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map) {
	clear();
	/* replace ids first (if needed) */
	if(ids_map) {
		unsigned int n = 0;
		for(it i = e; i != end; ++i) {
			auto x = ids_map->find(i->first);
			if(x == ids_map->end()) {
				ids_map->insert(std::make_pair(i->first,n));
				i->first = n;
				n++;
			}
			else i->first = x->second;
			x = ids_map->find(i->second);
			if(x == ids_map->end()) {
				ids_map->insert(std::make_pair(i->second,n));
				i->second = n;
				n++;
			}
			else i->second = x->second;
		}
	}
	std::sort(e,end,[](const std::pair<unsigned int,unsigned int> a, const std::pair<unsigned int,unsigned int> b) {
			return (a.first < b.first || (a.first == b.first && a.second < b.second) );
		});
	end = std::unique(e,end);
	return create_graph_sorted(e,end,0,true); /* we do not pass ids_map, it was already taken care of before */
}

/* 3. create a graph from edges supplied as iterators of std::pair<unsigned int,unsigned int>
 * 
 * the iterators should be forward iterators, ideally generating the data on the fly, as it will be
 * copied into temporary arrays */
template<class it, class sent>
int graph::create_graph_copy_sort(it& e, const sent& end, std::unordered_map<unsigned int,unsigned int>* ids_map) {
	clear();
	
	/* copy all data into the edges and temporary arrays */
	unsigned int* in_edges = 0;
	size_t in_edges_size = 0;
	size_t j = 0;
	int ret = 0;
	if( grow_edges() || grow_edges_s(&in_edges,&in_edges_size,EDGES_GROW) )
		{ fprintf(stderr,"graph::create_graph_copy_sort(): could not allocate memory!\n"); ret = 1; goto create_graph_copy_sort_err; }
	
	for(;e!=end;++e,j++) {
		if(j == edges_size) if(grow_edges()) { fprintf(stderr,"graph::create_graph_copy_sort(): could not allocate memory!\n"); ret = 1; goto create_graph_copy_sort_err; }
		if(j == in_edges_size) if(grow_edges_s(&in_edges,&in_edges_size,EDGES_GROW))
			{ fprintf(stderr,"graph::create_graph_copy_sort(): could not allocate memory!\n"); ret = 1; goto create_graph_copy_sort_err; }
		edges[j] = (*e).first;
		in_edges[j] = (*e).second;
	}
	nedges = j;
	ret = create_graph(in_edges, ids_map);
	
create_graph_copy_sort_err:
	if(in_edges) free(in_edges);
	return ret;
}		


#endif

