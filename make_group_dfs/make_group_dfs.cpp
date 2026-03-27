#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>


class Track_file {
public:
	int eventid, trackid, pl, rawid;
};
class Group_file :public Track_file {
public:
	int link_num;
	std::vector<std::tuple<int, int, int, int>> linklet;
};
class linklet_header {
public:
	int pos0, raw0, pos1, raw1;
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
size_t hash_value(const Segment& d)
{
	// 複数の値のハッシュ値を組み合わせてハッシュ値を計算するには、
	// boost::hash_combine を使います。
	size_t h = 0;
	boost::hash_combine(h, d.pos);
	boost::hash_combine(h, d.rawid);
	return h;
}

std::vector<Track_file> read_track(std::string filename);
void read_linklet_list(std::string filename, boost::unordered_multimap<Segment, Segment>& link);
Group_file make_group(Track_file& t, boost::unordered_multimap<Segment, Segment>& link);
void output_group(std::string filename, std::vector<Group_file>& g);
Group_file make_group_next(Track_file& t, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev);
Group_file make_group_prev(Track_file& t, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev);
Group_file group_merge(Track_file& t, Group_file& g_next, Group_file& g_prev);
void read_linklet_list(std::string filename, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev);
Group_file make_group_next_one(Track_file& t, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev);
bool judge_shower(Group_file& g, bool flg_multi);
int64_t Count_path(Group_file& g);


int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, "usage:prg file-in-track file_in_link file-out-track file-out-shower flg\n");
		fprintf(stderr, "flg=0:multi del\n");
		fprintf(stderr, "flg=1:multi \n");
		exit(1);
	}
	std::string file_in_track = argv[1];
	std::string file_in_link = argv[2];
	std::string file_out_gruop = argv[3];
	std::string file_out_gruop_shower = argv[4];
	bool flg_multi = std::stoi(argv[5]);

	std::vector<Track_file> track = read_track(file_in_track);
	boost::unordered_multimap<Segment, Segment> link_next;
	boost::unordered_multimap<Segment, Segment> link_prev;
	read_linklet_list(file_in_link, link_next, link_prev);

	std::vector<Group_file> group;
	std::vector<Group_file> group_shower;
	bool shower_flg = false;
	for (int i = 0; i < track.size(); i++) {
		fprintf(stderr, "make group %d/%d\r", i, track.size());
		Group_file g_next = make_group_next(track[i], link_next, link_prev);
		Group_file g_prev = make_group_prev(track[i], link_next, link_prev);
		Group_file merge = group_merge(track[i], g_next, g_prev);

		if (judge_shower(merge, flg_multi)) {
			//shower like pick up
			group_shower.push_back(merge);
			//printf("group%5d base %d link %d ratio %.3lf\n", merge.eventid, num_base, num_link, num_base*1. / num_link);
			//merge = make_group_next_one(track[i], link_next, link_prev);
		}
		else {
			group.push_back(merge);
		}

	}
	fprintf(stderr, "\n");

	output_group(file_out_gruop, group);
	if (group_shower.size() > 0) {
		output_group(file_out_gruop_shower, group_shower);
	}
	exit(0);
}
std::vector<Track_file> read_track(std::string filename) {
	std::ifstream ifs(filename);
	Track_file t;
	std::vector<Track_file> ret;
	int count = 0;
	while (ifs >> t.eventid >> t.trackid >> t.pl >> t.rawid) {
		if (count % 1000 == 0) {
			fprintf(stderr, "read track %d\r", count);
		}
		count++;
		ret.push_back(t);

	}
	fprintf(stderr, "read track %d\n", count);
	return ret;
}


int64_t Linklet_header_num(std::string filename) {
	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	return size2 / sizeof(linklet_header);
}

void read_linklet_list(std::string filename, boost::unordered_multimap<Segment, Segment>& link) {
	int64_t link_num = Linklet_header_num(filename);

	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	int64_t GB = size2 / (1000 * 1000 * 1000);
	int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
	int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
	if (GB > 0) {
		std::cout << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
	}
	else {
		std::cout << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
	}
	int64_t count = 0;
	linklet_header l;
	Segment seg0, seg1;

	while (ifs.read((char*)&l, sizeof(linklet_header))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;
		seg0.pos = l.pos0;
		seg0.rawid = l.raw0;
		seg1.pos = l.pos1;
		seg1.rawid = l.raw1;
		link.insert(std::make_pair(seg0, seg1));
		link.insert(std::make_pair(seg1, seg0));
	}
	auto size1 = eofpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
	if (count == 0) {
		fprintf(stderr, "%s no linklet!\n", filename.c_str());
		exit(1);
	}
}

void read_linklet_list(std::string filename, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev) {
	int64_t link_num = Linklet_header_num(filename);

	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	int64_t GB = size2 / (1000 * 1000 * 1000);
	int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
	int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
	if (GB > 0) {
		std::cerr << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
	}
	else {
		std::cerr << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
	}
	int64_t count = 0;
	linklet_header l;
	Segment seg0, seg1;

	while (ifs.read((char*)&l, sizeof(linklet_header))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;
		seg0.pos = l.pos0;
		seg0.rawid = l.raw0;
		seg1.pos = l.pos1;
		seg1.rawid = l.raw1;
		if (seg0.pos > seg1.pos)std::swap(seg0, seg1);

		link_next.insert(std::make_pair(seg0, seg1));
		link_prev.insert(std::make_pair(seg1, seg0));
	}
	auto size1 = eofpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
	if (count == 0) {
		fprintf(stderr, "%s no linklet!\n", filename.c_str());
		exit(1);
	}
}

Group_file make_group_next(Track_file& t, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev) {
	Group_file ret;
	ret.eventid = t.eventid;
	ret.trackid = t.trackid;
	ret.pl = t.pl;
	ret.rawid = t.rawid;

	Segment start, seg0, seg1;
	start.pos = t.pl * 10;
	start.rawid = t.rawid;



	std::set<std::pair<Segment, Segment>> all_link;
	std::set<Segment> seen;
	std::set<Segment> add, now;

	add.insert(start);
	int loop_num = 0;
	while (add.size() != 0) {
		now = add;
		add.clear();
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			seen.insert(*itr);
			//printf("loop=%d %d %d\n", loop_num, itr->pos, itr->rawid);
		}
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			if (link_next.find(*itr) == link_next.end())continue;
			auto range = link_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				if (seen.count(res->second) == 0) {
					add.insert(res->second);
				}
				if (res->first.pos < res->second.pos) {
					seg0.pos = res->first.pos;
					seg0.rawid = res->first.rawid;
					seg1.pos = res->second.pos;
					seg1.rawid = res->second.rawid;
				}
				else {
					seg1.pos = res->first.pos;
					seg1.rawid = res->first.rawid;
					seg0.pos = res->second.pos;
					seg0.rawid = res->second.rawid;
				}
				all_link.insert(std::make_pair(seg0, seg1));
			}
		}
		loop_num++;
	}

	seen.clear();
	for (auto itr = all_link.begin(); itr != all_link.end(); itr++) {
		add.insert(itr->first);
		add.insert(itr->second);
	}
	while (add.size() != 0) {
		now = add;
		add.clear();
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			seen.insert(*itr);
			//printf("loop=%d %d %d\n", loop_num, itr->pos, itr->rawid);
		}
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			if (link_prev.find(*itr) == link_prev.end())continue;
			auto range = link_prev.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				if (res->second.pos / 10 <= t.pl)continue;

				if (seen.count(res->second) == 0) {
					add.insert(res->second);
				}
				if (res->first.pos < res->second.pos) {
					seg0.pos = res->first.pos;
					seg0.rawid = res->first.rawid;
					seg1.pos = res->second.pos;
					seg1.rawid = res->second.rawid;
				}
				else {
					seg1.pos = res->first.pos;
					seg1.rawid = res->first.rawid;
					seg0.pos = res->second.pos;
					seg0.rawid = res->second.rawid;
				}
				all_link.insert(std::make_pair(seg0, seg1));
			}
		}
	}







	//printf("loop fin\n");
	//printf("link size=%d\n", all_link.size());
	if (all_link.size() > 0) {
		for (auto itr = all_link.begin(); itr != all_link.end(); itr++) {
			ret.linklet.push_back(std::make_tuple(itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid));
		}
	}
	ret.link_num = ret.linklet.size();
	return ret;
}

Group_file make_group_next_one(Track_file& t, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev) {
	Group_file ret;
	ret.eventid = t.eventid;
	ret.trackid = t.trackid;
	ret.pl = t.pl;
	ret.rawid = t.rawid;

	Segment start, seg0, seg1;
	start.pos = t.pl * 10;
	start.rawid = t.rawid;



	std::set<std::pair<Segment, Segment>> all_link;
	std::set<Segment> seen;
	std::set<Segment> add, now;

	add.insert(start);
	int loop_num = 0;
	while (add.size() != 0) {
		now = add;
		add.clear();
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			seen.insert(*itr);
			//printf("loop=%d %d %d\n", loop_num, itr->pos, itr->rawid);
		}
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			if (link_next.find(*itr) == link_next.end())continue;
			auto range = link_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				if (seen.count(res->second) == 0) {
					add.insert(res->second);
				}
				if (res->first.pos < res->second.pos) {
					seg0.pos = res->first.pos;
					seg0.rawid = res->first.rawid;
					seg1.pos = res->second.pos;
					seg1.rawid = res->second.rawid;
				}
				else {
					seg1.pos = res->first.pos;
					seg1.rawid = res->first.rawid;
					seg0.pos = res->second.pos;
					seg0.rawid = res->second.rawid;
				}
				all_link.insert(std::make_pair(seg0, seg1));
			}
		}
		loop_num++;
	}

	//printf("loop fin\n");
	//printf("link size=%d\n", all_link.size());
	if (all_link.size() > 0) {
		for (auto itr = all_link.begin(); itr != all_link.end(); itr++) {
			ret.linklet.push_back(std::make_tuple(itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid));
		}
	}
	ret.link_num = ret.linklet.size();
	return ret;
}

Group_file make_group_prev(Track_file& t, boost::unordered_multimap<Segment, Segment>& link_next, boost::unordered_multimap<Segment, Segment>& link_prev) {
	Group_file ret;
	ret.eventid = t.eventid;
	ret.trackid = t.trackid;
	ret.pl = t.pl;
	ret.rawid = t.rawid;

	Segment start, seg0, seg1;
	start.pos = t.pl * 10;
	start.rawid = t.rawid;



	std::set<std::pair<Segment, Segment>> all_link;
	std::set<Segment> seen;
	std::set<Segment> add, now;

	add.insert(start);
	while (add.size() != 0) {
		now = add;
		add.clear();
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			seen.insert(*itr);
			//printf("loop=%d %d %d\n", loop_num, itr->pos, itr->rawid);
		}
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			if (link_prev.find(*itr) == link_prev.end())continue;
			auto range = link_prev.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				if (seen.count(res->second) == 0) {
					add.insert(res->second);
				}
				if (res->first.pos < res->second.pos) {
					seg0.pos = res->first.pos;
					seg0.rawid = res->first.rawid;
					seg1.pos = res->second.pos;
					seg1.rawid = res->second.rawid;
				}
				else {
					seg1.pos = res->first.pos;
					seg1.rawid = res->first.rawid;
					seg0.pos = res->second.pos;
					seg0.rawid = res->second.rawid;
				}
				all_link.insert(std::make_pair(seg0, seg1));
			}
		}
	}

	seen.clear();
	for (auto itr = all_link.begin(); itr != all_link.end(); itr++) {
		add.insert(itr->first);
		add.insert(itr->second);
	}
	while (add.size() != 0) {
		now = add;
		add.clear();
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			seen.insert(*itr);
			//printf("loop=%d %d %d\n", loop_num, itr->pos, itr->rawid);
		}
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			if (link_next.find(*itr) == link_next.end())continue;
			auto range = link_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				if (res->second.pos / 10 >= t.pl)continue;

				if (seen.count(res->second) == 0) {
					add.insert(res->second);
				}
				if (res->first.pos < res->second.pos) {
					seg0.pos = res->first.pos;
					seg0.rawid = res->first.rawid;
					seg1.pos = res->second.pos;
					seg1.rawid = res->second.rawid;
				}
				else {
					seg1.pos = res->first.pos;
					seg1.rawid = res->first.rawid;
					seg0.pos = res->second.pos;
					seg0.rawid = res->second.rawid;
				}
				all_link.insert(std::make_pair(seg0, seg1));
			}
		}
	}







	//printf("loop fin\n");
	//printf("link size=%d\n", all_link.size());
	if (all_link.size() > 0) {
		for (auto itr = all_link.begin(); itr != all_link.end(); itr++) {
			ret.linklet.push_back(std::make_tuple(itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid));
		}
	}
	ret.link_num = ret.linklet.size();
	return ret;
}

Group_file group_merge(Track_file& t, Group_file& g_next, Group_file& g_prev) {
	Group_file ret;
	ret.eventid = t.eventid;
	ret.trackid = t.trackid;
	ret.pl = t.pl;
	ret.rawid = t.rawid;

	std::set<std::tuple<int, int, int, int>> linklet_all;
	for (int i = 0; i < g_next.linklet.size(); i++) {
		linklet_all.insert(g_next.linklet[i]);
	}
	for (int i = 0; i < g_prev.linklet.size(); i++) {
		linklet_all.insert(g_prev.linklet[i]);
	}

	if (linklet_all.size() == 0) {
		ret.link_num = 0;
		ret.linklet.clear();
		ret.linklet.shrink_to_fit();
		return ret;
	}
	for (auto itr = linklet_all.begin(); itr != linklet_all.end(); itr++) {
		ret.linklet.push_back(*itr);
	}
	ret.link_num = linklet_all.size();
	return ret;
}

void output_group(std::string filename, std::vector<Group_file>& g) {
	std::ofstream ofs(filename);
	int count = 0;
	for (auto itr = g.begin(); itr != g.end(); itr++) {
		if (count % 1000 == 0) {
			fprintf(stderr, "\r wrtie file %d/%d(%4.1lf%%)", count, g.size(), count * 100. / g.size());
		}
		count++;

		ofs << std::fixed << std::right
			<< std::setw(10) << std::setprecision(0) << itr->eventid << " "
			<< std::setw(5) << std::setprecision(0) << itr->trackid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pl << " "
			<< std::setw(4) << std::setprecision(0) << itr->rawid << " "
			<< std::setw(10) << std::setprecision(0) << itr->link_num << std::endl;
		if (itr->linklet.size() == 0)continue;
		for (auto itr2 = itr->linklet.begin(); itr2 != itr->linklet.end(); itr2++) {
			ofs << std::fixed << std::right
				<< std::setw(5) << std::setprecision(0) << std::get<0>(*itr2) << " "
				<< std::setw(5) << std::setprecision(0) << std::get<1>(*itr2) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<2>(*itr2) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<3>(*itr2) << std::endl;
		}
	}
	fprintf(stderr, "\r wrtie file %d/%d(%4.1lf%%)\n", count, g.size(), count * 100. / g.size());

}

bool judge_shower(Group_file& g, bool flg_multi) {

	std::set<std::pair<int, int>> basetrack;
	for (auto itr = g.linklet.begin(); itr != g.linklet.end(); itr++) {
		basetrack.insert(std::make_pair(std::get<0>(*itr), std::get<2>(*itr)));
		basetrack.insert(std::make_pair(std::get<1>(*itr), std::get<3>(*itr)));
	}
	int num_base = basetrack.size(), num_link = g.linklet.size();

	std::map<int, int> base_num;
	for (auto itr = basetrack.begin(); itr != basetrack.end(); itr++) {
		if (base_num.count(itr->first) == 0) {
			base_num.insert(std::make_pair(itr->first, 1));
		}
		else {
			base_num.at(itr->first) += 1;
		}
	}
	int max_base = 0;
	for (auto itr = base_num.begin(); itr != base_num.end(); itr++) {
		max_base = std::max(itr->second, max_base);
	}
	double rms, ratio;
	if (max_base == 0) {
		ratio = 1;
		rms = 0;
	}
	else {
		double sum = 0, sum2 = 0, count = 0;
		for (auto itr = base_num.begin(); itr != base_num.end(); itr++) {
			sum += itr->second;
			sum2 += itr->second * itr->second;
			count++;
		}
		if (count == 0) {
		}
		else {
			//1basetrackから平均何本linkletが出るか
			//1本のchainでは2に漸近
			ratio = num_link * 2. / num_base;
			if (sum2 / count - pow(sum / count, 2) < 0)rms = 0;
			else rms = sqrt(sum2 / count - pow(sum / count, 2));
		}
	}

	int64_t all_path_num = Count_path(g);


	bool ret = false;
	//multi消しなし
	if (flg_multi) {
		ret = (all_path_num > 1e8 || all_path_num < 0);
	}
	//multi消しあり
	else {
		ret = (all_path_num > 1e8 || all_path_num < 0);
	}

	printf("%5d %2d %5d %6d %4.3lf %4d %7.3lf %lld\n", g.eventid, ret, num_base, num_link, ratio, max_base, rms, all_path_num);

	return ret;

}

int64_t Count_path(Group_file& g) {
	boost::unordered_multimap<Segment, Segment> link_next;
	boost::unordered_multimap<Segment, Segment> link_prev;
	Segment seg0, seg1;
	std::map<Segment, int64_t*> seg_count;
	std::multimap<int, std::pair<Segment, int64_t>>seg_map;
	for (auto itr = g.linklet.begin(); itr != g.linklet.end(); itr++) {
		seg0.pos = std::get<0>(*itr);
		seg1.pos = std::get<1>(*itr);
		seg0.rawid = std::get<2>(*itr);
		seg1.rawid = std::get<3>(*itr);
		link_next.insert(std::make_pair(seg0, seg1));
		link_prev.insert(std::make_pair(seg1, seg0));
		seg_map.insert(std::make_pair(seg0.pos, std::make_pair(seg0, 0)));
		seg_map.insert(std::make_pair(seg1.pos, std::make_pair(seg1, 0)));
		//seg_map.insert(std::make_pair(seg1, 0));
	}
	for (auto itr = seg_map.begin(); itr != seg_map.end(); itr++) {
		if (link_prev.count(itr->second.first) == 0) {
			itr->second.second = 1;
		}
		seg_count.insert(std::make_pair(itr->second.first, &(itr->second.second)));
	}
	for (auto itr = seg_map.begin(); itr != seg_map.end(); itr++) {
		if (link_prev.count(itr->second.first) == 0) continue;
		auto range = link_prev.equal_range(itr->second.first);
		for (auto res = range.first; res != range.second; res++) {
			itr->second.second += *(seg_count.at(res->second));
		}
	}

	int64_t all_path_num = 0;
	for (auto itr = seg_map.begin(); itr != seg_map.end(); itr++) {
		if (link_next.count(itr->second.first) == 0) {
			if (INT64_MAX - all_path_num > itr->second.second) {
				all_path_num += itr->second.second;
			}
			else {
				all_path_num = -1;
				break;
			}
		}
	}
	return all_path_num;
}