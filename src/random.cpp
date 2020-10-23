/*
 * random.cpp
 *
 *  Created on: 22 Oct 2020
 *      Author: rich
 */

#include "random.h"
#include <experimental/random>
#include <unordered_set>

using namespace std;
using namespace tbb;

void RandomSubset(concurrent_vector<Interaction> & large, concurrent_vector<Interaction> & small){

	concurrent_vector<Interaction> replacement;
	int largeSize = large.size()-1;

	std::unordered_set<unsigned int> indices;
	while (indices.size() < small.size())
		indices.insert(std::experimental::fundamentals_v2::randint(0,largeSize));

	for (auto count:indices)
		replacement.push_back(large.at(count));

	large.swap(replacement);
}

void RandomChoose(concurrent_vector<Interaction> & first, concurrent_vector<Interaction> & second){
	if (first.size() > second.size())
	{
		if (second.size() < 0.9*first.size())
		{
			RandomSubset(first, second);
		}
	}
	else
	{
		if (first.size() < 0.9*second.size())
		{
			RandomSubset(second, first);
		}
	}
}
