#define NsegLimit 6


#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib, "lib_l2cx.lib")
#include <LibL2c-x.h>

#include <filesystem>
#include <set>
#include <map>
#include <netscan_data_types_ui.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
//#include <boost/unordered_set.hpp>
//#include <boost/unordered_map.hpp>
#include <stack>
#include <unordered_set>

class Track_file {
public:
	int eventid, trackid, pl, rawid;
};
class Group_file :public Track_file {
public:
	int link_num;
	std::vector<std::tuple<int, int, int, int>> linklet;
};
struct Segment {
	int16_t pos;
	int32_t rawid;
	bool operator==(const Segment& rhs) const
	{
		return pos == rhs.pos && rawid == rhs.rawid;
	}
	bool operator<(const Segment& rhs) const {
		if (pos == rhs.pos) {
			return rawid < rhs.rawid;
		}
		return pos < rhs.pos;
	}
};

namespace mfile0 {
	bool operator<(const M_Base& left, const M_Base& right) {
		if (left.pos == right.pos) {
			return left.rawid < right.rawid;
		}
		return left.pos < right.pos;
	}
	bool operator==(const M_Base& left, const M_Base& right) {
		return left.pos == right.pos && left.rawid == right.rawid;
	}
}
bool sort_chain(mfile0::M_Chain& left, mfile0::M_Chain& right) {
	if (left.basetracks.begin()->group_id != right.basetracks.begin()->group_id) {
		return left.basetracks.begin()->group_id < right.basetracks.begin()->group_id;
	}
	else if (left.basetracks.begin()->flg_i[1] != right.basetracks.begin()->flg_i[1]) {
		return left.basetracks.begin()->flg_i[1] < right.basetracks.begin()->flg_i[1];
	}
	else {
		return left.chain_id < right.chain_id;
	}
	return left.chain_id < right.chain_id;
}

std::vector <  Group_file > read_group_file(std::string filename);
std::map<std::pair<int, int>, std::vector<mfile0::M_Chain>> chain_divide_event(std::vector<mfile0::M_Chain>& chains);
mfile0::M_Chain select_chain(Group_file& g, std::vector<mfile0::M_Chain>& chains);
std::pair<Segment, Segment> single_path_connect(std::map<int, Segment>& chain_base, std::multimap<Segment, Segment>& path_next, std::multimap<Segment, Segment>& path_prev);
void Enumerate_Path(std::multimap<Segment, Segment>& path, std::vector<Segment >& hist, Segment start, int nseg_limit, std::vector<std::vector<Segment>>& all_path);
std::vector<std::vector<Segment>> make_nseg_chain(Segment& edge_seg, std::multimap<Segment, Segment>& path_next, int nseg_limit);
std::pair<double, Segment> judge_connect(Segment& edge_seg, std::map<int, Segment>& chain_base, std::multimap<Segment, Segment>& path_next, std::multimap<Segment, Segment>& path_prev, std::map<Segment, mfile0::M_Base>& base_map);
std::vector<Segment> target_chain_pickup(std::vector<std::vector<Segment>>& all_path, Segment& edge_seg, std::map<int, Segment>& chain_base);
double Calc_chain_value(std::vector<Segment>& seg_all, std::map<Segment, mfile0::M_Base>& base_map);
std::vector<std::vector<Segment>> multi_chain_to_unique(std::vector<std::vector<Segment>>& all_path);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "prg file-in-mfile file-in-group file-out-mfile\n");
		exit(1);
	}

	std::string file_in_mfile = argv[1];
	std::string file_in_group = argv[2];
	std::string file_out_mfile = argv[3];

	std::vector<Group_file>group = read_group_file(file_in_group);

	mfile0::Mfile m, out;
	mfile0::read_mfile(file_in_mfile, m);
	out.header = m.header;

	std::map<std::pair<int, int>, std::vector<mfile0::M_Chain>>chain_group = chain_divide_event(m.chains);
	std::pair<int, int>id;
	for (auto& g : group) {
		id.first = g.eventid;
		id.second = g.trackid;
		//printf( "groupid = %d chain id = %d\n", id.first, id.second);
		//printf("target segment %d,%d\n", g.pl, g.rawid);
		if (chain_group.count(id) == 0) {
			fprintf(stderr, "exception:groupid = %d chain id = %d\n", id.first, id.second);
		}

		auto chains = chain_group.at(id);
		//printf("chain num= %d\n", chains.size());

		mfile0::M_Chain chain_sel = select_chain(g, chains);
		out.chains.push_back(chain_sel);
	}

	std::sort(out.chains.begin(), out.chains.end(), sort_chain);
	mfile0::write_mfile(file_out_mfile, out);

}
std::vector <  Group_file > read_group_file(std::string filename) {
	std::vector <  Group_file > ret;
	std::ifstream ifs(filename);
	int gid, tid, muon_num, link_num, t_pl, t_rawid;
	int pl, raw;
	double weight;
	std::tuple<int, int, int, int> link;
	int count = 0;

	while (ifs >> gid >> tid >> t_pl >> t_rawid >> link_num) {
		printf("\r read group %d", count);
		count++;
		Group_file g;
		g.eventid = gid;
		g.trackid = tid;
		g.pl = t_pl;
		g.rawid = t_rawid;
		for (int i = 0; i < link_num; i++) {
			ifs >> std::get<0>(link) >> std::get<1>(link) >> std::get<2>(link) >> std::get<3>(link);
			g.linklet.push_back(link);
		}
		g.link_num = g.linklet.size();
		ret.push_back(g);
	}
	printf("\r read group %d\n", count);

	return ret;

}

std::map<std::pair<int, int>, std::vector<mfile0::M_Chain>> chain_divide_event(std::vector<mfile0::M_Chain>& chains) {
	std::multimap<std::pair<int, int>, mfile0::M_Chain> event_multimap;
	for (auto& c : chains) {
		event_multimap.insert(std::make_pair(std::make_pair(c.basetracks[0].group_id, c.basetracks[0].flg_i[1]), c));
	}

	std::map<std::pair<int, int>, std::vector<mfile0::M_Chain>> ret;
	for (auto itr = event_multimap.begin(); itr != event_multimap.end(); itr++) {

		std::vector<mfile0::M_Chain> c;
		auto range = event_multimap.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			c.push_back(res->second);
		}
		ret.insert(std::make_pair(itr->first, c));
		itr = std::next(itr, event_multimap.count(itr->first) - 1);
	}
	return ret;
}


mfile0::M_Chain select_chain(Group_file& g, std::vector<mfile0::M_Chain>& chains) {
	int pl = g.pl;
	int rawid = g.rawid;
	std::map<Segment, mfile0::M_Base> base_map;
	std::map<int, Segment> chain_base;

	std::set<std::pair<Segment, Segment>>path_set;
	std::multimap<Segment, Segment> path_next;
	std::multimap<Segment, Segment> path_prev;
	Segment seg0, seg1;
	//chainをbasetrackとlinkletに分解
	for (auto& c : chains) {
		for (int i = 0; i < c.basetracks.size(); i++) {
			seg1.pos = int(c.basetracks[i].pos / 10) * 10 + 1;
			seg1.rawid = c.basetracks[i].rawid;
			base_map.insert(std::make_pair(seg1, c.basetracks[i]));
			if (i != 0) {
				path_set.insert(std::make_pair(seg0, seg1));
			}
			seg0 = seg1;
		}
	}
	for (auto& path : path_set) {
		path_next.insert(std::make_pair(path.first, path.second));
		path_prev.insert(std::make_pair(path.second, path.first));
	}


	seg0.pos = g.pl * 10 + 1;
	seg0.rawid = g.rawid;
	chain_base.insert(std::make_pair(seg0.pos, seg0));
	int loop_num = 0;
	while (true) {
		//printf("loop %d\n", loop_num);
		loop_num += 1;
		std::pair<Segment, Segment> edge_pair = single_path_connect(chain_base, path_next, path_prev);
		//for (auto itr = chain_base.begin(); itr != chain_base.end(); itr++) {
		//	printf("PL =%3d rawid=%d\n", itr->second.pos / 10, itr->second.rawid);
		//}
		printf("edge 0\n");
		printf("PL =%3d rawid=%d\n", edge_pair.first.pos / 10, edge_pair.first.rawid);
		printf("edge 1\n");
		printf("PL =%3d rawid=%d\n", edge_pair.second.pos / 10, edge_pair.second.rawid);

		//上流
		std::pair<double, Segment> upstream_val = judge_connect(chain_base.rbegin()->second, chain_base, path_next, path_prev, base_map);
		//下流
		std::pair<double, Segment> downstream_val = judge_connect(chain_base.begin()->second, chain_base, path_prev, path_next, base_map);

		printf("up:%lf,down%lf\n", upstream_val.first, downstream_val.first);

		//上流下流両方つながる場合
		if (upstream_val.first > 0 && downstream_val.first > 0) {
			if (upstream_val.first > downstream_val.first) {
				chain_base.insert(std::make_pair(downstream_val.second.pos, downstream_val.second));
			}
			else {
				chain_base.insert(std::make_pair(upstream_val.second.pos, upstream_val.second));
			}
		}
		else if (upstream_val.first > 0) {
			chain_base.insert(std::make_pair(upstream_val.second.pos, upstream_val.second));
		}
		else if (downstream_val.first > 0) {
			chain_base.insert(std::make_pair(downstream_val.second.pos, downstream_val.second));
		}
		else {
			break;
		}
	}

	std::ofstream ofs("error_log_of_Result_multi_chain_partner.txt");
	ofs <<"PL\trawid" << std::endl;
	mfile0::M_Chain ret;
	for (auto& seg : chain_base) {
		auto res = base_map.find(seg.second);
		if (res == base_map.end()) {
			fprintf(stderr, "PL%03d rawid=%d not found\n", seg.second.pos / 10, seg.second.rawid);
			//exit(1);
			ofs << seg.second.pos / 10 << "\t" << seg.second.rawid << std::endl;
			std::cin.get();
			continue;
		}

		ret.basetracks.push_back(res->second);
	}
	ret.chain_id = g.trackid;
	ret.nseg = ret.basetracks.size();
	ret.pos0 = ret.basetracks.begin()->pos;
	ret.pos1 = ret.basetracks.rbegin()->pos;
	return ret;
}
//singleである限りpush back
std::pair<Segment, Segment> single_path_connect(std::map<int, Segment>& chain_base, std::multimap<Segment, Segment>& path_next, std::multimap<Segment, Segment>& path_prev) {
	std::pair<Segment, Segment> ret;
	Segment edge_seg;
	int next_count, prev_count;
	//上流-->下流
	while (true) {
		edge_seg = chain_base.rbegin()->second;
		next_count = path_next.count(edge_seg);
		if (next_count == 0) {
			ret.second.pos = -1;
			ret.second.rawid = -1;
			break;
		}
		if (next_count > 1) {
			ret.second = edge_seg;
			break;
		}
		auto next_seg = path_next.find(edge_seg);
		prev_count = path_prev.count(next_seg->second);
		if (prev_count > 1) {
			ret.second = edge_seg;
			break;
		}
		chain_base.insert(std::make_pair(next_seg->second.pos, next_seg->second));
	}
	//下流-->上流
	while (true) {
		edge_seg = chain_base.begin()->second;
		prev_count = path_prev.count(edge_seg);
		if (prev_count == 0) {
			ret.first.pos = -1;
			ret.first.rawid = -1;
			break;
		}
		if (prev_count > 1) {
			ret.first = edge_seg;
			break;
		}
		auto prev_seg = path_prev.find(edge_seg);
		next_count = path_next.count(prev_seg->second);
		if (next_count > 1) {
			ret.first = edge_seg;
			break;
		}
		chain_base.insert(std::make_pair(prev_seg->second.pos, prev_seg->second));
	}
	return ret;
}
//	
std::pair<double, Segment> judge_connect(Segment& edge_seg, std::map<int, Segment>& chain_base, std::multimap<Segment, Segment>& path_next, std::multimap<Segment, Segment>& path_prev, std::map<Segment, mfile0::M_Base>& base_map) {
	int nseg_limit = 6;

	std::pair<double, Segment> ret;
	ret.first = -1;
	ret.second.pos = -1;
	ret.second.rawid = -1;
	int count_next = path_next.count(edge_seg);

	if (count_next == 0)return ret;
	//逆側分岐(n-1)
	else if (count_next == 1) {
		auto next_seg = path_next.find(edge_seg);
		ret.second = next_seg->second;
		int count_prev = path_prev.count(next_seg->second);
		std::vector<std::vector<Segment>> all_path = make_nseg_chain(next_seg->second, path_prev, std::min(NsegLimit, int(chain_base.size() + 1)));

		all_path = multi_chain_to_unique(all_path);
		//edge segを含むものを分離
		//edge segを含み、chain baseに含まれないものを除去

		std::vector<Segment> target_chain = target_chain_pickup(all_path, edge_seg, chain_base);

		//seg-->basetrack情報の抽出
		//指標付け
		ret.first = Calc_chain_value(target_chain, base_map);
		//指標が大きい場合は他と比較
		//if (ret.first > 0.001) {
		for (auto& path : all_path) {
			if (ret.first > Calc_chain_value(path, base_map)) {
				ret.first = -1;
				break;
			}
		}

		//}
	}
	//1-n分岐
	else {
		auto range = path_next.equal_range(edge_seg);
		std::vector<std::vector<Segment>> all_path;
		std::vector<Segment> path_conv;
		if (chain_base.begin()->second == edge_seg) {
			for (auto itr = chain_base.begin(); itr != chain_base.end(); itr++) {
				path_conv.push_back(itr->second);
				if (path_conv.size() == NsegLimit) {
					break;
				}
			}
		}
		else if (chain_base.rbegin()->second == edge_seg) {
			for (auto itr = chain_base.rbegin(); itr != chain_base.rend(); itr++) {
				path_conv.push_back(itr->second);
				if (path_conv.size() == NsegLimit) {
					break;
				}
			}

		}
		else {
			fprintf(stderr, "exception direction not difine\n");
			exit(1);
		}

		for (auto res = range.first; res != range.second; res++) {
			std::vector<Segment> path_add = path_conv;
			path_add.push_back(res->second);
			all_path.push_back(path_add);
		}

		double val;
		for (auto& path : all_path) {
			val = Calc_chain_value(path, base_map);

			if (ret.first < 0 || val < ret.first) {
				ret.first = val;
				ret.second = *path.rbegin();
			}
		}
	}

	return ret;

}




//nseg分vector chainを生成
std::vector<std::vector<Segment>> make_nseg_chain(Segment& edge_seg, std::multimap<Segment, Segment>& path_next, int nseg_limit) {

	std::vector<std::vector<Segment>> all_path;
	std::vector<Segment >hist;
	Enumerate_Path(path_next, hist, edge_seg, nseg_limit, all_path);

	int nseg_min = nseg_limit;
	for (auto& path : all_path) {
		nseg_min = std::min(nseg_min, int(path.size()));
	}
	//printf("nseg min =%d\n", nseg_min);

	std::vector<std::vector<Segment>> ret;
	for (auto& path : all_path) {
		std::vector<Segment>one_path;

		for (auto& seg : path) {
			one_path.push_back(seg);
			if (one_path.size() == nseg_min) {
				ret.push_back(one_path);
				break;
			}
		}
	}

	for (int i = 0; i < ret.size(); i++) {
		printf("path id %d\n", i);
		for (int j = 0; j < ret[i].size(); j++) {
			printf("(%3d,%8d)--", ret[i][j].pos / 10, ret[i][j].rawid);
		}
		printf("\n");
	}

	return ret;
}
void Enumerate_Path(std::multimap<Segment, Segment>& path, std::vector<Segment >& hist, Segment start, int nseg_limit, std::vector<std::vector<Segment>>& all_path) {
	hist.push_back(start);
	//printf("[%d,%d]:add\n", start.pos,start.rawid);
	//for (int i = 0; i < hist.size(); i++) {
	//	if (i + 1 != hist.size()) {
	//		printf("(%d,%d)-->", hist[i].pos, hist[i].rawid);
	//	}
	//	else {
	//		printf("(%d,%d)\n", hist[i].pos, hist[i].rawid);
	//	}
	//}

	if (hist.size() - 1 == nseg_limit) {
		all_path.push_back(hist);
		hist.pop_back();
		return;
	}
	auto res = path.find(start);
	if (res == path.end()) {
		all_path.push_back(hist);
		hist.pop_back();
		return;
	}

	auto range = path.equal_range(start);
	for (auto itr = range.first; itr != range.second; itr++) {
		std::pair<Segment, Segment> next_path = std::make_pair(itr->first, itr->second);
		Enumerate_Path(path, hist, next_path.second, nseg_limit, all_path);
	}
	hist.pop_back();

}

std::vector<Segment> target_chain_pickup(std::vector<std::vector<Segment>>& all_path, Segment& edge_seg, std::map<int, Segment>& chain_base) {
	std::vector<Segment> ret;
	std::set<Segment> use_segment;
	for (auto itr = chain_base.begin(); itr != chain_base.end(); itr++) {
		use_segment.insert(itr->second);
	}
	int target_count = 0, remove_count = 0, remain_count = 0, add_count = 0;
	for (auto itr = all_path.begin(); itr != all_path.end();) {
		add_count = 0;

		bool flg = false;
		for (auto& seg : (*itr)) {
			if (use_segment.count(seg) == 1)flg = true;
			//if (edge_seg == seg)flg = true;
		}
		if (!flg) {
			remain_count += 1;
			itr++;
			continue;
		}

		for (auto& seg : (*itr)) {
			if (use_segment.count(seg) == 0) {
				add_count += 1;
			}
		}
		//printf("add count %d\n", add_count);
		//edgeの次のbasetrackのみaddしている場合は残す
		if (add_count == 1) {
			ret = *itr;
			target_count += 1;
		}
		itr = all_path.erase(itr);
		remove_count += 1;
	}
	//printf("remain path:%d\n", remain_count);
	//printf("target path:%d\n", target_count);
	//printf("remove path:%d\n", remove_count);

	return ret;

}
std::vector<std::vector<Segment>> multi_chain_to_unique(std::vector<std::vector<Segment>>& all_path) {
	std::vector<std::vector<Segment>>unique_path;
	for (auto& path : all_path) {
		bool flg = false;
		for (int i = 0; i < unique_path.size(); i++) {
			if (path == unique_path[i]) {
				flg = true;
			}
		}
		if (flg)continue;
		unique_path.push_back(path);
	}
	//printf("path unique %d --> %d\n", all_path.size(), unique_path.size());
	return unique_path;
}

//指標付け
double Calc_chain_value(std::vector<Segment>& seg_all, std::map<Segment, mfile0::M_Base>& base_map) {
	std::vector<mfile0::M_Base> basetracks;
	for (auto& seg : seg_all) {
		auto res = base_map.find(seg);
		if (res == base_map.end()) {
			fprintf(stderr, "PL%03d rawid=%d not found\n", seg.pos / 10, seg.rawid);
			exit(1);
		}
		basetracks.push_back(res->second);
	}
	//basetrackの位置から角度の計算

	//角度、角度ずれを計算
	double ax_chain, ay_chain;
	{
		int count = 0;
		double x = 0, y = 0, z = 0;
		double x2 = 0, y2 = 0, z2 = 0;
		double xz = 0, yz = 0;
		for (auto& b : basetracks) {
			count += 1;
			x += b.x;
			y += b.y;
			z += b.z;
			x2 += pow(b.x, 2);
			y2 += pow(b.y, 2);
			z2 += pow(b.z, 2);
			xz += b.x * b.z;
			yz += b.y * b.z;
			printf("ax:%.4lf, ay:%.4lf\n", b.ax, b.ay);

		}
		ax_chain = (count * xz - x * z) / (count * z2 - z * z);
		ay_chain = (count * yz - y * z) / (count * z2 - z * z);
	}
	printf("result\n");
	printf("ax:%.4lf, ay:%.4lf\n", ax_chain, ay_chain);

	//指標を出力
	double val = 0;
	for (auto& b : basetracks) {
		val += pow(b.ax - ax_chain, 2);
		val += pow(b.ay - ay_chain, 2);
	}
	//printf("value:%lf\n", val);
	return val;
}
