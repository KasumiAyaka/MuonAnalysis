#pragma once
#pragma once
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <fstream>
std::multimap<int, int> path_reverse2(std::multimap<int, int>&path);
void output_dot(std::string filename, std::set<int> &vertex, std::multimap<int, int> path_next);
std::set<int> Set_vertex(std::multimap<int, int>&path);
std::multimap<int, int>path_merge2(std::multimap<int, int>&path0, std::multimap<int, int>&path1);
void dfs_judge_path_exist(std::multimap<int, int>&path, std::vector<int>&hist, int start, int end, std::set<int>&seen, bool &flg);
void Enumerate_Path(std::multimap<int, int>&path, std::vector<int >&hist, int start, int goal, std::vector<std::vector<std::pair<int, int>>>&all_path);
void Enumerate_Cycle(std::multimap<int, int>&path, std::vector<std::vector<std::pair<int, int>>>&all_cycle);

std::vector<std::vector<std::pair<int, int>>> cycle_enumerate(std::multimap<int, int> path_next) {

	std::multimap<int, int> path_prev;
	std::multimap<int, int> path_all;
	std::set<int>vertex;

	path_prev = path_reverse2(path_next);
	//vertex = Set_vertex(path_next);
	path_all = path_merge2(path_next, path_prev);

	std::vector<std::vector<std::pair<int, int>>>all_cycle;
	Enumerate_Cycle(path_all, all_cycle);
	return all_cycle;


}
std::multimap<int, int> path_reverse2(std::multimap<int, int>&path) {
	std::multimap<int, int>  ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		ret.insert(std::make_pair(itr->second, itr->first));
	}
	return ret;
}
std::set<int> Set_vertex(std::multimap<int, int>&path) {
	std::set<int> ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		ret.insert(itr->first);
		ret.insert(itr->second);
	}
	return ret;

}
std::multimap<int, int>path_merge2(std::multimap<int, int>&path0, std::multimap<int, int>&path1) {
	std::multimap<int, int>ret;
	for (auto itr = path0.begin(); itr != path0.end(); itr++) {
		ret.insert(*itr);
	}
	for (auto itr = path1.begin(); itr != path1.end(); itr++) {
		ret.insert(*itr);
	}
	return ret;
}
void dfs_judge_path_exist(std::multimap<int, int>&path, std::vector<int>&hist, int start, int end, std::set<int>&seen, bool &flg) {
	seen.insert(start);
	hist.push_back(start);
	if (start == end) {
		flg = true;
		return;
	}
	if (path.count(start) == 0) {
		hist.pop_back();
		return;
	}
	auto range = path.equal_range(start);
	for (auto res = range.first; res != range.second; res++) {
		if (seen.count(res->second) == 1)continue;
		dfs_judge_path_exist(path, hist, res->second, end, seen, flg);
		if (flg)return;
	}
	hist.pop_back();
}
std::multimap<int, int> path_remove(std::multimap<int, int>&path, std::pair<int, int>&del_path) {
	std::multimap<int, int> ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		if (itr->first == del_path.first&&itr->second == del_path.second)continue;
		if (itr->first == del_path.second&&itr->second == del_path.first)continue;
		ret.insert(*itr);
	}
	return ret;
}
std::multimap<int, int> edge_remove(std::multimap<int, int>&path, int&del_edge) {
	std::multimap<int, int> ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		if (itr->first == del_edge || itr->second == del_edge)continue;
		ret.insert(*itr);
	}
	return ret;
}

void Enumerate_Path(std::multimap<int, int>&path, std::vector<int >&hist, int start, int goal, std::vector<std::vector<std::pair<int, int>>>&all_path) {
	hist.push_back(start);
	//printf("[%d]:add\n", start);
	//for (int i = 0; i < hist.size(); i++) {
	//	if (i + 1 != hist.size()) {
	//		printf("%d-->", hist[i]);
	//	}
	//	else {
	//		printf("%d\n", hist[i]);
	//	}
	//}

	if (start == goal) {
		//printf("detect path\n\t");
		//for (int i = 0; i < hist.size(); i++) {
		//	if (i + 1 != hist.size()) {
		//		printf("%d-->", hist[i]);
		//	}
		//	else {
		//		printf("%d\n", hist[i]);
		//	}
		//}
		std::vector<std::pair<int, int>> hist_recon;
		for (int i = 0; i < hist.size(); i++) {
			if (i + 1 != hist.size()) {
				hist_recon.push_back(std::make_pair(hist[i], hist[i + 1]));
			}
		}
		all_path.push_back(hist_recon);
		//hist.pop_back();
		return;
	}

	auto res = path.find(start);
	if (res == path.end()) {
		//hist.pop_back();
		return;
	}
	std::pair<int, int> next_path = std::make_pair(res->first, res->second);
	//printf("%d-%d:start\n", next_path.first, next_path.second);

	std::multimap<int, int>path2 = path_remove(path, next_path);
	bool flg = false;
	std::vector<int>hist2;
	std::set<int>seen;
	dfs_judge_path_exist(path2, hist2, start, goal, seen, flg);
	if (flg) {
		if (hist2.size() < 2) {
			fprintf(stderr, "path size small\n");
			fprintf(stderr, "function [ Enumerate_Path]\n");
			exit(1);
		}
		hist.pop_back();
		Enumerate_Path(path2, hist, start, goal, all_path);
	}
	//printf("%d-%d:end\n", next_path.first, next_path.second);
	path2 = edge_remove(path, start);
	Enumerate_Path(path2, hist, next_path.second, goal, all_path);
	hist.pop_back();
}

void Enumerate_Cycle(std::multimap<int, int>&path, std::vector<std::vector<std::pair<int, int>>>&all_cycle) {

	while (path.size() > 0) {
		std::pair<int, int> cut_path = *path.begin();
		path = path_remove(path, cut_path);
		std::vector<int >hist2;
		std::vector<std::vector<std::pair<int, int>>>all_path;
		//next_pathÇêÿÇ¡ÇΩÇ§Ç¶Ç≈pathÇ™Ç†ÇÈÇ©
		std::set<int>seen;
		bool flg = false;
		dfs_judge_path_exist(path, hist2, cut_path.first, cut_path.second, seen, flg);
		if (flg) {
			hist2.clear();
			Enumerate_Path(path, hist2, cut_path.first, cut_path.second, all_path);
			//printf("cut:%d -%d\n", cut_path.first, cut_path.second);
			for (int i = 0; i < all_path.size(); i++) {
				std::vector<std::pair<int, int>>cycle;
				cycle.push_back(std::make_pair(cut_path.second, cut_path.first));
				for (int j = 0; j < all_path[i].size(); j++) {
					cycle.push_back(all_path[i][j]);
				}
				all_cycle.push_back(cycle);
			}
		}
		Enumerate_Cycle(path, all_cycle);
	}
}

void output_dot(std::string filename, std::set<int> &vertex, std::multimap<int, int> path_next) {
	std::ofstream ofs(filename);
	ofs << "graph{" << std::endl;
	for (auto itr = vertex.begin(); itr != vertex.end(); itr++) {
		if (std::next(itr, 1) == vertex.end()) {
			ofs << *itr << std::endl;
		}
		else {
			ofs << *itr << ",";
		}
	}
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		ofs << itr->first << "--" << itr->second << std::endl;
	}
	ofs << "}" << std::endl;

}

