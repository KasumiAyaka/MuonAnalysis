// mfile-->txt
// kasumi
// 2024/08/11
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Btrk {
	int pl;
	int64_t rawid;
};
bool operator<(const Btrk& lhs, const Btrk& rhs) {
	return std::tie(lhs.pl, lhs.rawid) < std::tie(rhs.pl, rhs.rawid);
}


std::vector<mfile0::M_Chain>track_selection_PL(std::vector<mfile0::M_Chain>& chain, int pl, int num_pl_up, int num_pl_down);
std::vector<mfile0::M_Chain>chain_dlat_selection(std::vector<mfile0::M_Chain>& chain, double threshold);
void output_base_rawid(std::string filename, std::vector<mfile0::M_Chain>& chain, int pl);
std::multimap<Btrk, int> list(std::string filename) ;
std::vector<mfile0::M_Chain>Search_basetrack(std::vector<mfile0::M_Chain> & chain, std::multimap<Btrk, int> & list);
void output(std::string filename, std::vector<mfile0::M_Chain>& chain);

int main(int argc, char** argv) {

	if (argc != 4) {
		fprintf(stderr, "exception\n");
		fprintf(stderr, "usage:in-mfile.bmf list.txt output.txt\n");
		exit(1);
	}

	std::string filename_mfile = argv[1];
	std::string filename_list = argv[2];
	std::string filename_output = argv[3];

	mfile0::Mfile m;

	mfile1::read_mfile_extension(filename_mfile, m);

	std::multimap<Btrk, int> l=list(filename_list);
	std::vector<mfile0::M_Chain> chain=Search_basetrack(m.chains, l);
	output(filename_output, chain);


}
std::vector<mfile0::M_Chain>track_selection_PL(std::vector<mfile0::M_Chain>& chain, int pl, int num_pl_up, int num_pl_down) {

	std::vector<mfile0::M_Chain> ret;


	int hit_up_num, hit_down_num, base_pl, pl_distance;
	int64_t count = 0;
	for (int64_t i = 0; i < chain.size(); i++) {
		if (count % 10000 == 0) {
			printf("\r track selection PL%03d %lld/%lld(%4.1lf%%)", pl, count, chain.size(), count * 100. / chain.size());
		}
		count++;


		hit_down_num = 0;
		hit_up_num = 0;
		for (int j = 0; j < chain[i].basetracks.size(); j++) {
			base_pl = chain[i].basetracks[j].pos / 10;
			pl_distance = abs(base_pl - pl);
			//downstream
			if (base_pl < pl) {
				if (pl_distance <= num_pl_down)hit_down_num++;
			}
			//upstream
			else if (base_pl > pl) {
				if (pl_distance <= num_pl_up)hit_up_num++;
			}

		}
		if (hit_up_num == num_pl_up && hit_down_num == num_pl_down) {
			ret.push_back(chain[i]);
		}


	}
	printf("\r track selection PL%03d %lld/%lld(%4.1lf%%)\n", pl, count, chain.size(), count * 100. / chain.size());

	printf("chain num %lld -->%lld (%4.1lf%%)\n", chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;
}
std::vector<mfile0::M_Chain>chain_dlat_selection(std::vector<mfile0::M_Chain>& chain, double threshold) {
	std::vector<mfile0::M_Chain> ret;
	double ax, ay;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		ax = mfile0::chain_ax(*itr);
		ay = mfile0::chain_ay(*itr);
		//â^ìÆó ï™ïzÇÇ›ÇÈ
		//printf("",mfile0::angle_diff_dev_lat(*itr, ax, ay));
		if (mfile0::angle_diff_dev_lat(*itr, ax, ay) > threshold)continue;
		ret.push_back(*itr);
	}
	printf("chain lateral selection <= %5.4lf : %d--> %d (%4.1lf%%)\n", threshold, chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;

}

void output_base_rawid(std::string filename, std::vector<mfile0::M_Chain>& chain, int mode) {

	//std::ofstream ofs(filename);
	//int base_pl;
	//int64_t count = 0;
	//for (int64_t i = 0; i < chain.size(); i++) {
	//	if (count % 10000 == 0) {
	//		printf("\r output base rawid %lld/%lld(%4.1lf%%)", count, chain.size(), count * 100. / chain.size());
	//	}
	//	count++;
	//	for (int j = 0; j < chain[i].basetracks.size(); j++) {
	//		base_pl = chain[i].basetracks[j].pos / 10;
	//		if (base_pl == pl) {
	//			ofs << std::right << std::fixed
	//				<< std::setw(5) << base_pl << ""
	//				<< std::setw(12) << chain[i].basetracks[j].rawid << std::endl;
	//		}
	//	}
	//}

	std::ofstream ofs(filename);
	if (mode == 0) {
		std::cout << "output: flg chainid pl groupid rawid ph ax ay x y z" << std::endl;
		int count = 0;
		for (auto itr = chain.begin(); itr != chain.end(); itr++) {
			//Ç¬Ç»Ç™Ç¡ÇΩbasetrackï™ïz
			// chainç≈è„ó¨Å`ç≈â∫ó¨àÍå¬éËëOÇ‹Ç≈ÇÕflg=1,ç≈â∫ó¨ÇÕflg=0ÇèoóÕÇ∑ÇÈ-------------------------
			int flg = 0;
			for (int k = 0; k < itr->basetracks.size(); k++) {
				decltype(itr->basetracks)::iterator itr3 = std::next(itr->basetracks.begin(), k);
				// std::cout << itr3->ax << std::endl;

				flg = 0;
				if (k == itr->basetracks.size() - 1) {
					flg = 1;
					ofs << std::right << std::fixed
						<< std::setw(3) << std::setprecision(0) << flg << " "
						<< std::setw(20) << std::setprecision(0) << itr->chain_id << " "
						<< std::setw(4) << std::setprecision(0) << itr3->pos / 10 << " "
						<< std::setw(20) << std::setprecision(0) << itr3->group_id << " "
						<< std::setw(12) << std::setprecision(0) << itr3->rawid << " "
						<< std::setw(6) << std::setprecision(0) << itr3->ph << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ax << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ay << " "
						<< std::setw(8) << std::setprecision(1) << itr3->x << " "
						<< std::setw(8) << std::setprecision(1) << itr3->y << " "
						<< std::setw(8) << std::setprecision(1) << itr3->z << std::endl;
				}
				else {
					flg = 0;
					ofs << std::right << std::fixed
						<< std::setw(3) << std::setprecision(0) << flg << " "
						<< std::setw(20) << std::setprecision(0) << itr->chain_id << " "
						<< std::setw(4) << std::setprecision(0) << itr3->pos / 10 << " "
						<< std::setw(20) << std::setprecision(0) << itr3->group_id << " "
						<< std::setw(12) << std::setprecision(0) << itr3->rawid << " "
						<< std::setw(6) << std::setprecision(0) << itr3->ph << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ax << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ay << " "
						<< std::setw(8) << std::setprecision(1) << itr3->x << " "
						<< std::setw(8) << std::setprecision(1) << itr3->y << " "
						<< std::setw(8) << std::setprecision(1) << itr3->z << std::endl;
				}


			}
			count++;
			//if (count == 10)break;
			//--- chainç≈è„ó¨Å`ç≈â∫ó¨àÍå¬éËëOÇ‹Ç≈ÇÕflg=1,ç≈â∫ó¨ÇÕflg=0ÇèoóÕÇ∑ÇÈ--- Ç±Ç±Ç‹Ç≈------------
		}

	}
	if (mode == 1) {
		std::cout << "output: chainid pl0 pl1 groupid rawid ph phv ax ay x y z count nseg" << std::endl;
		for (auto itr = chain.begin(); itr != chain.end(); itr++) {
			// average of ph,ax,ay
			// ---out put : C-id pl G-id R-id ave-ph ave-ax ave-ay x,y,z-------------------------
			int flg = 0;
			int count = 0;
			double bph = 0.;
			double bax = 0.;
			double bay = 0.;
			double bvph = 0.;
			for (int k = 0; k < itr->basetracks.size(); k++) {// chainÇç\ê¨ÇµÇƒÇ¢ÇÈtrackÇ…Ç¬Ç¢ÇƒëñÇÈ
				auto itr0 = std::next(itr->basetracks.begin(), k);
				if (k == itr->basetracks.size() - 1) {
					count++;
					bph = bph + itr0->ph / 10000;
					bax = bax + itr0->ax;
					bay = bay + itr0->ay;
					bvph = bvph + itr0->ph % 10000;

					bph = bph / count;
					bax = bax / count;
					bay = bay / count;
					bvph = bvph / count;

					ofs << std::right << std::fixed
						<< std::setw(20) << std::setprecision(0) << itr->chain_id << " "
						<< std::setw(4) << std::setprecision(0) << itr->pos0 / 10 << " "
						<< std::setw(4) << std::setprecision(0) << itr->pos1 / 10 << " "
						<< std::setw(20) << std::setprecision(0) << itr0->group_id << " "
						<< std::setw(12) << std::setprecision(0) << itr0->rawid << " "
						<< std::setw(6) << std::setprecision(0) << bph << " "
						<< std::setw(6) << std::setprecision(0) << bvph << " "
						<< std::setw(7) << std::setprecision(4) << bax << " "
						<< std::setw(7) << std::setprecision(4) << bay << " "
						<< std::setw(8) << std::setprecision(1) << itr0->x << " " // most downstream
						<< std::setw(8) << std::setprecision(1) << itr0->y << " "
						<< std::setw(8) << std::setprecision(1) << itr0->z << " "
						<< count << " " << itr->nseg << std::endl;
				}
				else {
					count++;
					bph = bph + itr0->ph / 10000;
					bax = bax + itr0->ax;
					bay = bay + itr0->ay;
					bvph = bvph + itr0->ph % 10000;
				}
			}
		}
	}
	else exit(1);
	//int mod = 1;
	//if (mod == 0) {
	//	//int64_t count = 0;
	////	std::ofstream ofs(filename);
	//	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
	//		ofs << std::right << std::fixed
	//			<< std::setw(20) << std::setprecision(0) << itr->chain_id << " "
	//			<< std::setw(4) << std::setprecision(0) << itr->nseg << " "
	//			<< std::setw(4) << std::setprecision(0) << itr->pos0 / 10 << " "
	//			<< std::setw(4) << std::setprecision(0) << itr->pos1 / 10 << " "
	//			<< std::setw(20) << std::setprecision(0) << itr->basetracks.begin()->group_id << " "
	//			<< std::setw(12) << std::setprecision(0) << itr->basetracks.begin()->rawid << " "
	//			<< std::setw(6) << std::setprecision(0) << itr->basetracks.begin()->ph << " "
	//			<< std::setw(7) << std::setprecision(4) << itr->basetracks.begin()->ax << " "
	//			<< std::setw(7) << std::setprecision(4) << itr->basetracks.begin()->ay << " "
	//			<< std::setw(8) << std::setprecision(1) << itr->basetracks.begin()->x << " "
	//			<< std::setw(8) << std::setprecision(1) << itr->basetracks.begin()->y << " "
	//			<< std::setw(8) << std::setprecision(1) << itr->basetracks.begin()->z << std::endl;
	//	}

	//	//printf("\r output base rawid %lld/%lld(%4.1lf%%)\n", count, chain.size(), count * 100. / chain.size());
	//}
}

std::multimap<Btrk, int> list(std::string filename) {
	std::multimap<Btrk, int> ret;

	std::ifstream ifs;
	if (!ifs) {
		std::cerr << "File open error\nFilename : " << filename << std::endl;
		exit(1);
	}

	Btrk b;
	while (ifs >> b.pl >> b.rawid) {
		ret.insert(std::make_pair(b, b.rawid));
	}

	return ret;
}
std::vector<mfile0::M_Chain>Search_basetrack(std::vector<mfile0::M_Chain>& chain, std::multimap<Btrk, int>& list){

	std::vector<mfile0::M_Chain> ret;
	Btrk b;

	int hit_up_num, hit_down_num, base_pl, pl_distance;
	int64_t count = 0;
	for (int64_t i = 0; i < chain.size(); i++) {
		if (count % 10000 == 0) {
			printf("\r track search %lld/%lld(%4.1lf%%)",count, chain.size(), count * 100. / chain.size());
		}
		count++;

		for (int j = 0; j < chain[i].basetracks.size(); j++) {

			b.pl = chain[i].basetracks[j].pos / 10;
			b.rawid = chain[i].basetracks[j].rawid;

			if (list.find(b) != list.end()) {
				std::cout << "Find : " <<std::setw(3)<< b.pl << " " <<std::setw(20)<< b.rawid << std::endl;
				ret.push_back(chain[i]);

			}
		}


	}
	printf("\r track search %lld/%lld(%4.1lf%%)\n", count, chain.size(), count * 100. / chain.size());
	printf("chain num %lld -->%lld (%4.1lf%%)\n", chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;
}
void output(std::string filename, std::vector<mfile0::M_Chain>& chain) {


	std::ofstream ofs(filename);
		for (auto itr = chain.begin(); itr != chain.end(); itr++) {
			//Ç¬Ç»Ç™Ç¡ÇΩbasetrackï™ïz
			// chainç≈è„ó¨Å`ç≈â∫ó¨àÍå¬éËëOÇ‹Ç≈ÇÕflg=1,ç≈â∫ó¨ÇÕflg=0ÇèoóÕÇ∑ÇÈ-------------------------
			int flg = 0;
			for (int k = 0; k < itr->basetracks.size(); k++) {
				decltype(itr->basetracks)::iterator itr3 = std::next(itr->basetracks.begin(), k);
				// std::cout << itr3->ax << std::endl;

				flg = 0;
				if (k == itr->basetracks.size() - 1) {
					flg = 1;
					ofs << std::right << std::fixed
						//<< std::setw(3) << std::setprecision(0) << flg << " "
						<< std::setw(20) << std::setprecision(0) << itr->chain_id << " "
						<< std::setw(4) << std::setprecision(0) << itr3->pos / 10 << " "
						<< std::setw(20) << std::setprecision(0) << itr3->group_id << " "
						<< std::setw(12) << std::setprecision(0) << itr3->rawid << " "
						<< std::setw(6) << std::setprecision(0) << itr3->ph << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ax << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ay << " "
						<< std::setw(8) << std::setprecision(1) << itr3->x << " "
						<< std::setw(8) << std::setprecision(1) << itr3->y << " "
						<< std::setw(8) << std::setprecision(1) << itr3->z << std::endl;
				}
				else {
					flg = 0;
					ofs << std::right << std::fixed
						//<< std::setw(3) << std::setprecision(0) << flg << " "
						<< std::setw(20) << std::setprecision(0) << itr->chain_id << " "
						<< std::setw(4) << std::setprecision(0) << itr3->pos / 10 << " "
						<< std::setw(20) << std::setprecision(0) << itr3->group_id << " "
						<< std::setw(12) << std::setprecision(0) << itr3->rawid << " "
						<< std::setw(6) << std::setprecision(0) << itr3->ph << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ax << " "
						<< std::setw(7) << std::setprecision(4) << itr3->ay << " "
						<< std::setw(8) << std::setprecision(1) << itr3->x << " "
						<< std::setw(8) << std::setprecision(1) << itr3->y << " "
						<< std::setw(8) << std::setprecision(1) << itr3->z << std::endl;
				}


			}
			//count++;
			//if (count == 10)break;
			//--- chainç≈è„ó¨Å`ç≈â∫ó¨àÍå¬éËëOÇ‹Ç≈ÇÕflg=1,ç≈â∫ó¨ÇÕflg=0ÇèoóÕÇ∑ÇÈ--- Ç±Ç±Ç‹Ç≈------------
		}

}


namespace mfile0 {
	void mfile0::read_mfile(std::string file_path, Mfile& mfile) {
		std::ifstream ifs(file_path);
		//filesizeéÊìæ
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
			std::cout << "FILE size :" << GB << "." << std::setw(3) << std::setfill('0') << MB << " [GB]" << std::endl;
		}
		else {
			std::cout << "FILE size :" << MB << "." << std::setw(3) << std::setfill('0') << KB << " [MB]" << std::endl;
		}

		std::string str;

		M_Header head_tmp;
		for (int i = 0; i < 3; i++) {
			std::getline(ifs, str);
			head_tmp.head[i] = str;
		}
		std::getline(ifs, str);
		head_tmp.num_all_plate = stoi(str);

		{
			std::getline(ifs, str);
			auto str_v = StringSplit(str);
			for (int i = 0; i < str_v.size(); i++) {
				head_tmp.all_pos.push_back(stoi(str_v[i]));
			}
		}

		mfile.header = head_tmp;
		int cnt = 0;

		while (std::getline(ifs, str)) {
			if (cnt % 10000 == 0) {
				nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			cnt++;
			M_Chain chains;

			auto strs = StringSplit(str);
			if (strs.size() == 4) {
				chains.chain_id = stoi(strs[0]);
				chains.nseg = stoi(strs[1]);
				chains.pos0 = stoi(strs[2]);
				chains.pos1 = stoi(strs[3]);

				for (int i = 0; i < chains.nseg; i++) {

					M_Base bt;

					std::getline(ifs, str);
					strs = StringSplit(str);
					if (strs.size() == 15) {
						bt.pos = stoi(strs[0]);
						bt.group_id = stoi(strs[1]);
						bt.rawid = stoi(strs[2]);
						bt.ph = stoi(strs[3]);
						bt.ax = stof(strs[4]);
						bt.ay = stof(strs[5]);
						bt.x = stod(strs[6]);
						bt.y = stod(strs[7]);
						bt.z = stod(strs[8]);

						bt.flg_i[0] = stoi(strs[9]);
						bt.flg_i[1] = stoi(strs[10]);
						bt.flg_i[2] = stoi(strs[11]);
						bt.flg_i[3] = stoi(strs[12]);
						bt.flg_d[0] = stod(strs[13]);
						bt.flg_d[1] = stod(strs[14]);
					}
					else {
						fprintf(stderr, "mfile format err! cannot read file\n");
						fprintf(stderr, "base track information 15 -->%d\n", int(strs.size()));
						for (int i = 0; i < strs.size(); i++) {
							printf("%s ", strs[i].c_str());
						}
						printf("\n");
						throw std::exception();
					}
					chains.basetracks.push_back(bt);
				}
			}
			else {
				fprintf(stderr, "mfile format err! cannot read file\n");
				fprintf(stderr, "chain information 4 -->%d\n", int(strs.size()));
				for (int i = 0; i < strs.size(); i++) {
					printf("%s ", strs[i].c_str());
				}
				printf("\n");

				throw std::exception();
			}
			mfile.chains.push_back(chains);
		}
		auto size1 = eofpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
		if (cnt == 0) {
			fprintf(stderr, "%s no chain!\n", file_path.c_str());
			exit(1);
		}
	}
	void mfile0::read_mfile(std::string file_path, Mfile& mfile, int nseg_thr) {
		std::ifstream ifs(file_path);
		//filesizeéÊìæ
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
			std::cout << "FILE size :" << GB << "." << std::setw(3) << std::setfill('0') << MB << " [GB]" << std::endl;
		}
		else {
			std::cout << "FILE size :" << MB << "." << std::setw(3) << std::setfill('0') << KB << " [MB]" << std::endl;
		}

		std::string str;

		M_Header head_tmp;
		for (int i = 0; i < 3; i++) {
			std::getline(ifs, str);
			head_tmp.head[i] = str;
		}
		std::getline(ifs, str);
		head_tmp.num_all_plate = stoi(str);

		{
			std::getline(ifs, str);
			auto str_v = StringSplit(str);
			for (int i = 0; i < str_v.size(); i++) {
				head_tmp.all_pos.push_back(stoi(str_v[i]));
			}
		}

		mfile.header = head_tmp;
		int cnt = 0;

		while (std::getline(ifs, str)) {
			if (cnt % 10000 == 0) {
				nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			cnt++;
			M_Chain chains;

			auto strs = StringSplit(str);
			if (strs.size() == 4) {
				chains.chain_id = stoi(strs[0]);
				chains.nseg = stoi(strs[1]);
				chains.pos0 = stoi(strs[2]);
				chains.pos1 = stoi(strs[3]);
				if (chains.nseg < nseg_thr) {
					for (int i = 0; i < chains.nseg; i++) {
						std::getline(ifs, str);
						continue;
					}
				}
				else {
					for (int i = 0; i < chains.nseg; i++) {

						M_Base bt;

						std::getline(ifs, str);
						strs = StringSplit(str);
						if (strs.size() == 15) {
							bt.pos = stoi(strs[0]);
							bt.group_id = stoi(strs[1]);
							bt.rawid = stoi(strs[2]);
							bt.ph = stoi(strs[3]);
							bt.ax = stof(strs[4]);
							bt.ay = stof(strs[5]);
							bt.x = stod(strs[6]);
							bt.y = stod(strs[7]);
							bt.z = stod(strs[8]);

							bt.flg_i[0] = stoi(strs[9]);
							bt.flg_i[1] = stoi(strs[10]);
							bt.flg_i[2] = stoi(strs[11]);
							bt.flg_i[3] = stoi(strs[12]);
							bt.flg_d[0] = stod(strs[13]);
							bt.flg_d[1] = stod(strs[14]);
						}
						else {
							fprintf(stderr, "mfile format err! cannot read file\n");
							fprintf(stderr, "base track information 15 -->%d\n", int(strs.size()));
							for (int i = 0; i < strs.size(); i++) {
								printf("%s ", strs[i].c_str());
							}
							printf("\n");
							throw std::exception();
						}
						chains.basetracks.push_back(bt);
					}
				}
			}
			else {
				fprintf(stderr, "mfile format err! cannot read file\n");
				fprintf(stderr, "chain information 4 -->%d\n", int(strs.size()));
				for (int i = 0; i < strs.size(); i++) {
					printf("%s ", strs[i].c_str());
				}
				printf("\n");

				throw std::exception();
			}
			mfile.chains.push_back(chains);
		}
		auto size1 = eofpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
		if (cnt == 0) {
			fprintf(stderr, "%s no chain!\n", file_path.c_str());
			exit(1);
		}
	}

	void mfile0::write_mfile(std::string filename, const Mfile& mfile, int output) {
		std::ofstream ofs(filename);
		if (!ofs) {
			//file open s
			fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
			throw std::exception();
		}
		write_mfile_header(ofs, mfile.header, output);
		int count = 0;
		for (auto itr = mfile.chains.begin(); itr != mfile.chains.end(); itr++) {
			if (output == 1 && count % 10000 == 0) {
				fprintf(stderr, "\r Write Mfile chain ... %d/%d (%4.1lf%%)", count, int(mfile.chains.size()), count * 100. / mfile.chains.size());
			}
			count++;
			write_mfile_chain(ofs, *itr);
		}
		if (output == 1) {
			fprintf(stderr, "\r Write Mfile chain ... %d/%d (%4.1lf%%)\n", count, int(mfile.chains.size()), count * 100. / mfile.chains.size());
		}
	}
	void mfile0::write_mfile_header(std::ofstream& ofs, const M_Header& header, int output) {
		std::cout << std::right << std::fixed;
		ofs << header.head[0] << std::endl;
		ofs << header.head[1] << std::endl;
		ofs << header.head[2] << std::endl;
		ofs << std::setw(6) << header.num_all_plate << std::endl;
		for (int i = 0; i < header.all_pos.size(); i++) {
			ofs << " " << std::setw(4) << header.all_pos[i];
		}
		ofs << std::endl;
		if (output == 1) {
			fprintf(stderr, "mfile header write fin\n");
		}
	}
	void mfile0::write_mfile_chain(std::ofstream& ofs, const M_Chain& chains) {
		ofs << std::right << std::fixed
			<< std::setw(10) << chains.chain_id << " "
			<< std::setw(5) << chains.nseg << " "
			<< std::setw(6) << chains.pos0 << " "
			<< std::setw(6) << chains.pos1 << std::endl;
		for (auto& p : chains.basetracks) {
			ofs << std::right << std::fixed
				<< std::setw(6) << p.pos << " "
				<< std::setw(20) << p.group_id << " "
				<< std::setw(9) << p.rawid << " "
				<< std::setw(9) << p.ph << " "
				<< std::setw(8) << std::setprecision(4) << p.ax << " "
				<< std::setw(8) << std::setprecision(4) << p.ay << " "
				<< std::setw(10) << std::setprecision(1) << p.x << " "
				<< std::setw(10) << std::setprecision(1) << p.y << " "
				<< std::setw(10) << std::setprecision(0) << p.z << " "
				<< std::setw(2) << p.flg_i[0] << " "
				<< std::setw(2) << p.flg_i[1] << " "
				<< std::setw(2) << p.flg_i[2] << " "
				<< std::setw(2) << p.flg_i[3] << " "
				<< std::setw(5) << std::setprecision(4) << p.flg_d[0] << " "
				<< std::setw(5) << std::setprecision(4) << p.flg_d[1] << " "
				<< std::endl;
		}
	}
	double mfile0::chain_ax(mfile0::M_Chain chain) {
		double tmp = 0;
		int count = 0;
		for (auto itr = chain.basetracks.begin(); itr != chain.basetracks.end(); itr++) {
			tmp += itr->ax;
			count++;
		}
		return tmp / count;
	}
	double mfile0::chain_ay(mfile0::M_Chain chain) {
		double tmp = 0;
		int count = 0;
		for (auto itr = chain.basetracks.begin(); itr != chain.basetracks.end(); itr++) {
			tmp += itr->ay;
			count++;
		}
		return tmp / count;
	}
	double mfile0::angle_diff_dev_rad(M_Chain chain, double ax, double ay) {
		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double angle = sqrt(ax * ax + ay * ay);
		int count = 0;
		for (auto itr = chain.basetracks.begin(); itr != chain.basetracks.end(); itr++) {
			diff = (itr->ax * ax + itr->ay * ay) / angle - angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
		}
		return sqrt(sum2 / count - pow(sum / count, 2.));
	}
	double mfile0::angle_diff_dev_lat(M_Chain chain, double ax, double ay) {
		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double angle = sqrt(ax * ax + ay * ay);
		int count = 0;
		for (auto itr = chain.basetracks.begin(); itr != chain.basetracks.end(); itr++) {
			diff = (itr->ay * ax - itr->ax * ay) / angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
		}
		return sqrt(sum2 / count - pow(sum / count, 2.));
	}
	void apply_vph_correction(mfile0::Mfile& m, std::string filename_corr) {
		struct Vol_Correction {
			std::map<int, double> mean, sigma;
		};
		Vol_Correction ang_00_02, ang_02_03, ang_03_06, ang_06_08, ang_08_10, ang_10_12, ang_12_;

		std::ifstream ifs(filename_corr);
		std::string str;

		while (std::getline(ifs, str)) {
			auto str_v = StringSplit(str);
			if (str_v.size() == 15) {
				int pl = stoi(str_v[0]);
				ang_00_02.mean.insert(std::make_pair(pl, stod(str_v[1])));
				ang_02_03.mean.insert(std::make_pair(pl, stod(str_v[3])));
				ang_03_06.mean.insert(std::make_pair(pl, stod(str_v[5])));
				ang_06_08.mean.insert(std::make_pair(pl, stod(str_v[7])));
				ang_08_10.mean.insert(std::make_pair(pl, stod(str_v[9])));
				ang_10_12.mean.insert(std::make_pair(pl, stod(str_v[11])));
				ang_12_.mean.insert(std::make_pair(pl, stod(str_v[13])));

				ang_00_02.sigma.insert(std::make_pair(pl, stod(str_v[2])));
				ang_02_03.sigma.insert(std::make_pair(pl, stod(str_v[4])));
				ang_03_06.sigma.insert(std::make_pair(pl, stod(str_v[6])));
				ang_06_08.sigma.insert(std::make_pair(pl, stod(str_v[8])));
				ang_08_10.sigma.insert(std::make_pair(pl, stod(str_v[10])));
				ang_10_12.sigma.insert(std::make_pair(pl, stod(str_v[12])));
				ang_12_.sigma.insert(std::make_pair(pl, stod(str_v[14])));
			}
			else {
				fprintf(stderr, "incorect file format [%s]\n", filename_corr.c_str());
				exit(1);
			}

		}

		int pl_tmp;
		double angle_tmp;
		int count = 0;
		for (int i = 0; i < m.chains.size(); i++) {
			if (count % 100000 == 0) {
				fprintf(stderr, "\r vph correction apply %d/%d(%4.1lf%%)", count, int(m.chains.size()), count * 100. / m.chains.size());
			}
			count++;

			for (auto itr = m.chains[i].basetracks.begin(); itr != m.chains[i].basetracks.end(); itr++) {

				pl_tmp = itr->pos / 10;
				angle_tmp = sqrt(pow(itr->ax, 2.) + pow(itr->ay, 2.));
				if (angle_tmp < 0.2) itr->flg_d[0] = (itr->ph % 10000) / ang_00_02.mean[pl_tmp];
				else if (0.2 <= angle_tmp && angle_tmp < 0.3) itr->flg_d[0] = (itr->ph % 10000) / ang_02_03.mean[pl_tmp];
				else if (0.3 <= angle_tmp && angle_tmp < 0.6) itr->flg_d[0] = (itr->ph % 10000) / ang_03_06.mean[pl_tmp];
				else if (0.6 <= angle_tmp && angle_tmp < 0.8) itr->flg_d[0] = (itr->ph % 10000) / ang_06_08.mean[pl_tmp];
				else if (0.8 <= angle_tmp && angle_tmp < 1.0) itr->flg_d[0] = (itr->ph % 10000) / ang_08_10.mean[pl_tmp];
				else if (1.0 <= angle_tmp && angle_tmp < 1.2) itr->flg_d[0] = (itr->ph % 10000) / ang_10_12.mean[pl_tmp];
				else itr->flg_d[0] = (itr->ph % 10000) / ang_12_.mean[pl_tmp];
			}
		}

		fprintf(stderr, "\r vph correction apply %d/%d(%4.1lf%%)\n", count, int(m.chains.size()), count * 100. / m.chains.size());

	}
	void write_chaininf(std::string filename, std::vector<M_Chain_inf> chain_inf) {
		std::ofstream ofs(filename);

		if (!ofs) {
			//file open é∏îs
			fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
			exit(1);
		}
		if (chain_inf.size() == 0) {
			fprintf(stderr, "target Chain Information ... null\n");
			fprintf(stderr, "File[%s] has no text\n", filename.c_str());
		}
		else {
			int count = 0;
			for (auto itr = chain_inf.begin(); itr != chain_inf.end(); itr++) {
				if (count % 10000 == 0) {
					fprintf(stderr, "\r Write chain information dump ... %d/%d (%4.1lf%%)", count, int(chain_inf.size()), count * 100. / chain_inf.size());
				}
				count++;
				ofs << std::right << std::fixed
					<< std::setw(10) << std::setprecision(0) << itr->chainID << " "
					<< std::setw(10) << std::setprecision(0) << itr->groupID << " "
					<< std::setw(4) << std::setprecision(0) << itr->nseg << " "
					<< std::setw(4) << std::setprecision(0) << itr->pos0 << " "
					<< std::setw(4) << std::setprecision(0) << itr->pos1 << " "
					<< std::setw(7) << std::setprecision(4) << itr->ax << " "
					<< std::setw(7) << std::setprecision(4) << itr->ay << " "

					<< std::setw(10) << std::setprecision(1) << itr->x_down.second << " "
					<< std::setw(10) << std::setprecision(1) << itr->y_down.second << " "
					<< std::setw(10) << std::setprecision(1) << itr->x_down.first << " "
					<< std::setw(2) << std::setprecision(0) << itr->outflg_down << " "

					<< std::setw(10) << std::setprecision(1) << itr->x_up.second << " "
					<< std::setw(10) << std::setprecision(1) << itr->y_up.second << " "
					<< std::setw(10) << std::setprecision(1) << itr->x_up.first << " "
					<< std::setw(2) << std::setprecision(0) << itr->outflg_up << " "

					<< std::setw(7) << std::setprecision(4) << itr->vph_ratio << " "
					<< std::setw(7) << std::setprecision(4) << itr->vph_slope << " "
					<< std::setw(7) << std::setprecision(4) << itr->vph_slope_acc << " "

					<< std::setw(6) << std::setprecision(4) << itr->radial_deviation << " "
					<< std::setw(6) << std::setprecision(4) << itr->lateral_deviation << " "
					<< std::setw(6) << std::setprecision(4) << itr->d_radial_deviation << " "
					<< std::setw(6) << std::setprecision(4) << itr->d_lateral_deviation << std::endl;
			}
			fprintf(stderr, "\r Write chain information dump ... %d/%d (%4.1lf%%)\n", count, int(chain_inf.size()), count * 100. / chain_inf.size());

		}

		ofs.close();
	}
	M_Chain_inf chain2inf(mfile0::M_Chain c) {
		M_Chain_inf ret;
		ret.ax = mfile0::chain_ax(c);
		ret.ay = mfile0::chain_ay(c);
		ret.chainID = c.chain_id;
		ret.groupID = c.basetracks[0].group_id;
		ret.nseg = c.nseg;
		ret.pos0 = c.pos0;
		ret.pos1 = c.pos1;
		ret.x_down = std::make_pair(c.basetracks[0].z, c.basetracks[0].x);
		ret.y_down = std::make_pair(c.basetracks[0].z, c.basetracks[0].y);
		ret.x_up = std::make_pair(c.basetracks[c.nseg - 1].z, c.basetracks[c.nseg - 1].x);
		ret.y_up = std::make_pair(c.basetracks[c.nseg - 1].z, c.basetracks[c.nseg - 1].y);
		ret.outflg_down = false;
		ret.outflg_up = false;

		double  buf[4];
		std::pair<int, double> sum[4], sum2[4];
		for (int i = 0; i < 4; i++) {
			sum[i].first = 0;
			sum[i].second = 0;
			sum2[i].first = 0;
			sum2[i].second = 0;
		}
		buf[0] = ret.ax;
		buf[1] = ret.ay;
		double diff[2];
		//äpìxè¨Ç≥Ç¢Ç∆Ç´ÇÕx-y
		if (sqrt(pow(ret.ax, 2.) + pow(ret.ay, 2.)) < 0.01) {
			for (int i = 0; i < c.basetracks.size(); i++) {
				diff[0] = buf[0] - c.basetracks[i].ax;
				diff[1] = buf[1] - c.basetracks[i].ay;
				sum[0].first++;
				sum2[0].first++;
				sum[0].second += diff[0];
				sum2[0].second += pow(diff[0], 2.);
				sum[1].first++;
				sum2[1].first++;
				sum[1].second += diff[1];
				sum2[1].second += pow(diff[1], 2.);
				if (i != 0) {
					diff[0] = buf[2] - c.basetracks[i].ax;
					diff[1] = buf[3] - c.basetracks[i].ay;

					sum[2].first++;
					sum2[2].first++;
					sum[2].second += diff[0];
					sum2[2].second += pow(diff[0], 2.);
					sum[3].first++;
					sum2[3].first++;
					sum[3].second += diff[1];
					sum2[3].second += pow(diff[1], 2.);
				}
				buf[2] = c.basetracks[i].ax;
				buf[3] = c.basetracks[i].ay;
			}

		}
		//ëÂÇ´Ç¢Ç∆Ç´ÇÕrad-lat
		else {
			for (int i = 0; i < c.basetracks.size(); i++) {
				{
					double ax_axis = buf[0];
					double ay_axis = buf[1];
					double ax = c.basetracks[i].ax;
					double ay = c.basetracks[i].ay;

					double r = sqrt(pow(ax_axis, 2) + pow(ay_axis, 2));

					diff[0] = r - (ax_axis * ax + ay_axis * ay) / r;
					diff[1] = (ax_axis * ay - ay_axis * ax) / r;

				}
				//radial äpìxç∑
				sum[0].first++;
				sum2[0].first++;
				sum[0].second += diff[0];
				sum2[0].second += pow(diff[0], 2.);
				//lateral äpìxç∑
				sum[1].first++;
				sum2[1].first++;
				sum[1].second += diff[1];
				sum2[1].second += pow(diff[1], 2.);
				if (i != 0) {
					{
						double ax_axis = buf[2];
						double ay_axis = buf[3];
						double ax = c.basetracks[i].ax;
						double ay = c.basetracks[i].ay;

						double r = sqrt(pow(ax_axis, 2) + pow(ay_axis, 2));

						diff[0] = r - (ax_axis * ax + ay_axis * ay) / r;
						diff[1] = (ax_axis * ay - ay_axis * ax) / r;

					}
					//radial äpìxç∑
					sum[2].first++;
					sum2[2].first++;
					sum[2].second += diff[0];
					sum2[2].second += pow(diff[0], 2.);
					//lateral äpìxç∑
					sum[3].first++;
					sum2[3].first++;
					sum[3].second += diff[1];
					sum2[3].second += pow(diff[1], 2.);
				}
				buf[2] = c.basetracks[i].ax;
				buf[3] = c.basetracks[i].ay;
			}
		}
		double ave, ave_2;
		for (int i = 0; i < 4; i++) {
			if (sum[i].first != 0) {
				ave = sum[i].second / sum[i].first;
				ave_2 = sum2[i].second / sum2[i].first;
				if (ave_2 - pow(ave, 2) >= 0) {
					switch (i)
					{
					case 0:
						ret.radial_deviation = sqrt(ave_2 - pow(ave, 2));
						break;
					case 1:
						ret.lateral_deviation = sqrt(ave_2 - pow(ave, 2));
						break;
					case 2:
						ret.d_radial_deviation = sqrt(ave_2 - pow(ave, 2));
						break;
					case 3:
						ret.d_lateral_deviation = sqrt(ave_2 - pow(ave, 2));
						break;
					default:
						break;
					}
				}
				//ä€ÇﬂåÎç∑Ç≈É}ÉCÉiÉXÇ…
				else if (ave_2 - pow(ave, 2) > -1 * pow(0.001, 2)) {
					switch (i)
					{
					case 0:
						ret.radial_deviation = 0;
						break;
					case 1:
						ret.lateral_deviation = 0;
						break;
					case 2:
						ret.d_radial_deviation = 0;
						break;
					case 3:
						ret.d_lateral_deviation = 0;
						break;
					default:
						break;
					}
				}
				else {
					fprintf(stderr, "Chain ID %lld  cannot calc(deviation)\n", ret.chainID);
					fprintf(stderr, "i=%d in root [%lf]\n", i, ave_2 - pow(ave, 2));
					exit(1);
				}
			}
			else {
				switch (i)
				{
				case 0:
					ret.radial_deviation = -1;
					break;
				case 1:
					ret.lateral_deviation = -1;
					break;
				case 2:
					ret.d_radial_deviation = -1;
					break;
				case 3:
					ret.d_lateral_deviation = -1;
					break;
				default:
					break;
				}

			}

		}

		//vphÇÃåvéZ
		//ç≈è¨ìÒèÊñ@
		double x, xx, xy, y, num;
		double slope, slope_err, intercept, intercept_err, mean;
		x = 0;
		xx = 0;
		xy = 0;
		y = 0;
		num = 0;
		for (auto itr = c.basetracks.begin(); itr != c.basetracks.end(); itr++) {
			x += double(itr->pos / 10);
			xx += pow(double(itr->pos / 10), 2.);
			xy += double(itr->pos / 10) * itr->flg_d[0];
			y += itr->flg_d[0];
			num++;
		}
		if (num <= 1) {
			intercept = 0;
			intercept_err = 0;
			slope = 0;
			slope_err = 0;
			mean = 0;
		}
		else {
			double delta = num * xx - pow(x, 2);
			intercept = (xx * y - x * xy) / delta;
			slope = (num * xy - x * y) / delta;
			mean = y / num;
			//åÎç∑ÇÃå©êœÇ‡ÇË
			//ìùåv2à»â∫Ç≈ÇÕíºê¸ÇÕàÍà”Ç…íËÇ‹ÇÈÇΩÇﬂåvéZïsâ¬
			if (num <= 2) {
				intercept_err = 0;
				slope_err = 0;
			}
			else {
				double err_y = 0;
				for (auto itr = c.basetracks.begin(); itr != c.basetracks.end(); itr++) {
					err_y += pow(itr->flg_d[0] - (intercept + slope * double(itr->pos / 10)), 2);
				}
				err_y = sqrt(err_y / (num - 2));
				intercept_err = err_y * sqrt(xx / delta);
				slope_err = err_y * sqrt(num / delta);
			}
		}
		ret.vph_slope = slope;
		ret.vph_slope_acc = slope_err;
		ret.vph_ratio = mean;

		return ret;
	}
	void chain_inf_flg(M_Chain_inf& inf, mfile0::M_Chain c, int pl0, int pl1, double range[4], std::map<int, double> z) {
		if (c.pos0 / 10 <= pl0) {
			inf.outflg_down = true;
		}
		else {
			double position[2];
			int pl = c.basetracks[0].pos / 10;
			position[0] = c.basetracks[0].x + c.basetracks[0].ax * (z[pl - 1] - z[pl]);
			position[1] = c.basetracks[0].y + c.basetracks[0].ay * (z[pl - 1] - z[pl]);
			if (position[0] < range[0] || range[1] < position[0] || position[1] < range[2] || range[3] < position[1]) {
				inf.outflg_down = true;
			}
		}
		if (c.pos1 / 10 >= pl1) {
			inf.outflg_up = true;
		}
		else {
			double position[2];
			int pl = c.basetracks[0].pos / 10;
			position[0] = c.basetracks[c.nseg - 1].x + c.basetracks[c.nseg - 1].ax * (z[pl + 1] - z[pl]);
			position[1] = c.basetracks[c.nseg - 1].y + c.basetracks[c.nseg - 1].ay * (z[pl + 1] - z[pl]);
			if (position[0] < range[0] || range[1] < position[0] || position[1] < range[2] || range[3] < position[1]) {
				inf.outflg_up = true;
			}
		}

	}
	void set_header(int pl0, int pl1, mfile0::Mfile& m) {
		mfile0::M_Header head;
		for (int pl = pl0; pl <= pl1; pl++) {
			head.all_pos.push_back(pl * 10 + 1);
		}
		head.num_all_plate = int(head.all_pos.size());
		head.head[0] = "% Created by mkmf";
		head.head[1] = "% on 2007 / 12 / 25 23:28 : 09 + 09:00 (JST)";
		head.head[2] = "0       0   3   0      0.0   0.0000";
		m.header = head;
	}
	void mfile0::write_mfile_chain_IVE(std::ofstream& ofs, const M_Chain& chains) {
		//IVEï`âÊópèoóÕ
		double base_extra = 180;
		for (int i = 0; i < chains.basetracks.size(); i++) {
			ofs << std::right << std::fixed
				<< std::setw(12) << chains.chain_id << " "
				<< std::setw(10) << chains.basetracks.begin()->group_id << " "
				<< std::setw(4) << chains.basetracks[i].pos / 10 << " "
				<< std::setw(12) << chains.basetracks[i].rawid << " "
				<< std::setw(2) << chains.basetracks[i].ph / 10000 << " "
				<< std::setw(4) << chains.basetracks[i].ph % 10000 << " "
				//basetrackï`âÊparam
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].x - chains.basetracks[i].ax * base_extra << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].y - chains.basetracks[i].ay * base_extra << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].z - base_extra << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].x + chains.basetracks[i].ax * base_extra << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].y + chains.basetracks[i].ay * base_extra << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].z + base_extra << " "
				//chainï`âÊparam
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].x << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].y << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[i].z << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[std::max(0, i - 1)].x << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[std::max(0, i - 1)].y << " "
				<< std::setw(8) << std::setprecision(1) << chains.basetracks[std::max(0, i - 1)].z << std::endl;
		}
	}
	double chain_fit_dist(M_Chain chain, double slope[3], double intercept[3], double& dist) {
		if (chain.nseg < 3) {
			fprintf(stderr, "segment < 3 cannot fit\n");
			exit(1);
		}
		std::vector<double> x_v, y_v, z_v;
		for (auto itr = chain.basetracks.begin(); itr != chain.basetracks.end(); itr++) {
			x_v.push_back(itr->x);
			y_v.push_back(itr->y);
			z_v.push_back(itr->z);
		}
		double slope_error[2] = {};
		double position_error[2] = {};
		//z-xÇ≈ÇÃç≈è¨ìÒèÊ
		double x, y, xx, xy, n;
		x = 0;
		y = 0;
		xx = 0;
		xy = 0;
		n = 0;
		for (int i = 0; i < z_v.size(); i++) {
			x += z_v[i];
			y += x_v[i];
			xx += z_v[i] * z_v[i];
			xy += z_v[i] * x_v[i];
			n++;
		}
		slope[0] = (n * xy - x * y) / (n * xx - x * x);
		intercept[0] = (xx * y - xy * x) / (n * xx - x * x);
		//fitÇ©ÇÁÇÃç∑
		for (int i = 0; i < z_v.size(); i++) {
			slope_error[0] += pow((slope[0] * z_v[i] + intercept[0]) - x_v[i], 2);
		}
		position_error[0] = sqrt(slope_error[0] / (n - 2));
		slope_error[0] = position_error[0] * sqrt(n / (n * xx - x * x));

		//z-yÇ≈ÇÃç≈è¨ìÒèÊ
		x = 0;
		y = 0;
		xx = 0;
		xy = 0;
		n = 0;
		for (int i = 0; i < z_v.size(); i++) {
			x += z_v[i];
			y += y_v[i];
			xx += z_v[i] * z_v[i];
			xy += z_v[i] * y_v[i];
			n++;
		}
		slope[1] = (n * xy - x * y) / (n * xx - x * x);
		intercept[1] = (xx * y - xy * x) / (n * xx - x * x);
		//fitÇ©ÇÁÇÃç∑
		for (int i = 0; i < z_v.size(); i++) {
			slope_error[1] += pow((slope[1] * z_v[i] + intercept[1]) - y_v[i], 2);
		}
		position_error[1] = sqrt(slope_error[1] / (n - 2));
		slope_error[1] = position_error[1] * sqrt(n / (n * xx - x * x));

		slope[2] = 1;
		intercept[2] = 0;

		double sigma = 3;
		//à íuÇÃÇŒÇÁÇ¬Ç´
		dist = sqrt(pow(position_error[0] * sigma, 2) + pow(position_error[1] * sigma, 2));
		//åXÇ´Ç…ëŒÇ∑ÇÈåÎç∑Ç™è¨Ç≥Ç¢Ç‡ÇÃÇëIÇ‘
		return sqrt(pow(slope_error[0], 2) + pow(slope_error[1], 2));
	}
}

namespace mfile1 {
	//mfile txt-->binary
	void mfile1::converter(const mfile0::Mfile& old, MFile& mfile) {
		int count = 0;
		for (auto& c : old.chains) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r mfile convert txt-->binay (%4.1lf%%)", count * 100. / old.chains.size());
			}
			count++;

			MFileChain1 chain;
			chain.chain_id = c.chain_id;
			chain.nseg = c.nseg;
			chain.pos0 = c.pos0;
			chain.pos1 = c.pos1;
			mfile.chains.emplace_back(chain);
			std::vector< MFileBase1> basetracks;
			for (auto& b : c.basetracks) {
				MFileBase1 base;

				base.pos = b.pos;
				base.group_id = b.group_id;
				base.rawid = b.rawid;
				base.ph = b.ph;
				base.ax = b.ax;
				base.ay = b.ay;
				base.x = b.x;
				base.y = b.y;
				base.z = b.z;
				basetracks.emplace_back(base);
			}
			mfile.all_basetracks.emplace_back(basetracks);
		}
		fprintf(stderr, "\r mfile convert txt-->binay (%4.1lf%%)\n", count * 100. / old.chains.size());
	}
	//mfile binaty-->txt
	void mfile1::converter(const MFile& old, mfile0::Mfile& mfile) {
		mfile.header.head[0] = "% Created by mkmf";
		mfile.header.head[1] = "% on 2007 / 12 / 25 23:28 : 09 + 09:00 (JST)";
		mfile.header.head[2] = "	0       0   3   0      0.0   0.0000";
		std::set<int> pos;

		int count = 0;
		for (int i = 0; i < old.chains.size(); i++) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r mfile convert binay-->txt (%4.1lf%%)", count * 100. / old.chains.size());
			}
			count++;

			mfile0::M_Chain chain;
			chain.chain_id = old.chains[i].chain_id;
			chain.nseg = old.chains[i].nseg;
			chain.pos0 = old.chains[i].pos0;
			chain.pos1 = old.chains[i].pos1;
			for (int j = 0; j < old.all_basetracks[i].size(); j++) {
				mfile0::M_Base base;
				base.pos = old.all_basetracks[i][j].pos;
				base.group_id = old.all_basetracks[i][j].group_id;
				base.rawid = old.all_basetracks[i][j].rawid;
				base.ph = old.all_basetracks[i][j].ph;
				base.ax = old.all_basetracks[i][j].ax;
				base.ay = old.all_basetracks[i][j].ay;
				base.x = old.all_basetracks[i][j].x;
				base.y = old.all_basetracks[i][j].y;
				base.z = old.all_basetracks[i][j].z;
				base.flg_d[0] = 0;
				base.flg_d[1] = 0;
				base.flg_i[0] = 0;
				base.flg_i[1] = 0;
				base.flg_i[2] = 0;
				base.flg_i[3] = 0;

				pos.insert(base.pos);
				chain.basetracks.push_back(base);
			}
			mfile.chains.push_back(chain);
		}
		fprintf(stderr, "\r mfile convert binay-->txt (%4.1lf%%)\n", count * 100. / old.chains.size());

		mfile.header.num_all_plate = int(pos.size());

		for (auto itr = pos.begin(); itr != pos.end(); itr++) {
			mfile.header.all_pos.push_back(*itr);
		}
	}
	//mfile binaty-->txt
	void mfile1::converter(const MFile_minimum& old, mfile0::Mfile& mfile) {
		mfile.header.head[0] = "% Created by mkmf";
		mfile.header.head[1] = "% on 2007 / 12 / 25 23:28 : 09 + 09:00 (JST)";
		mfile.header.head[2] = "	0       0   3   0      0.0   0.0000";
		std::set<int> pos;

		int count = 0;
		for (int i = 0; i < old.chains.size(); i++) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r mfile convert binay-->txt (%4.1lf%%)", count * 100. / old.chains.size());
			}
			count++;

			mfile0::M_Chain chain;
			chain.chain_id = old.chains[i].chain_id;
			chain.nseg = old.chains[i].nseg;
			chain.pos0 = old.chains[i].pos0;
			chain.pos1 = old.chains[i].pos1;
			for (int j = 0; j < old.all_basetracks[i].size(); j++) {
				mfile0::M_Base base;
				base.pos = old.all_basetracks[i][j].pos;
				base.group_id = old.all_basetracks[i][j].group_id;
				base.rawid = old.all_basetracks[i][j].rawid;
				base.ph = old.all_basetracks[i][j].ph;
				base.ax = old.all_basetracks[i][j].ax;
				base.ay = old.all_basetracks[i][j].ay;
				base.x = old.all_basetracks[i][j].x;
				base.y = old.all_basetracks[i][j].y;
				base.z = old.all_basetracks[i][j].z;
				base.flg_d[0] = 0;
				base.flg_d[1] = 0;
				base.flg_i[0] = 0;
				base.flg_i[1] = 0;
				base.flg_i[2] = 0;
				base.flg_i[3] = 0;

				pos.insert(base.pos);
				chain.basetracks.push_back(base);
			}
			mfile.chains.push_back(chain);
		}
		fprintf(stderr, "\r mfile convert binay-->txt (%4.1lf%%)\n", count * 100. / old.chains.size());

		mfile.header.num_all_plate = int(pos.size());

		for (auto itr = pos.begin(); itr != pos.end(); itr++) {
			mfile.header.all_pos.push_back(*itr);
		}
	}

	void mfile1::read_mfile(std::string filepath, MFile& mfile) {
		std::ifstream ifs(filepath, std::ios::binary);
		//filesizeéÊìæ
		ifs.seekg(0, std::ios::end);
		int64_t eofpos = ifs.tellg();
		ifs.clear();
		ifs.seekg(0, std::ios::beg);
		int64_t begpos = ifs.tellg();
		int64_t size2 = eofpos - begpos;
		int64_t GB = size2 / (1000 * 1000 * 1000);
		int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
		int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
		if (GB > 0) {
			std::cout << "FILE size :" << GB << "." << std::setw(3) << std::setfill('0') << MB << " [GB]" << std::endl;
		}
		else {
			std::cout << "FILE size :" << MB << "." << std::setw(3) << std::setfill('0') << KB << " [MB]" << std::endl;
		}
		//Mfile headerÇÃì«Ç›çûÇ›
		ifs.read((char*)&mfile.header, sizeof(MFileHeader));
		if (ifs.eof()) { throw std::exception(); }
		std::string  filetype = "mfile-a0";
		memcpy((char*)filetype.data(), (char*)&mfile.header.filetype, filetype.size());
		if (filetype != "mfile-a0") { throw std::exception("File format is not mfile-a0."); }

		//mfile info headerÇÃì«Ç›çûÇ›
		ifs.read((char*)&mfile.info_header, sizeof(MFileInfoHeader));
		if (ifs.eof()) { throw std::exception(); }

		if (sizeof(MFileChain) != mfile.info_header.classsize1) { throw std::exception("Classsize1 is wrong."); }
		if (sizeof(MFileBase) != mfile.info_header.classsize2) { throw std::exception("Classsize2 is wrong."); }

		std::vector< MFileChain1> chains;
		chains.reserve(mfile.info_header.Nchain);

		std::vector< std::vector< MFileBase1>> all_basetracks;
		uint64_t count = 0;
		for (uint64_t c = 0; c < mfile.info_header.Nchain; c++) {
			if (count % 100000 == 0) {
				auto nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			count++;
			MFileChain1 chain;
			ifs.read((char*)&chain, sizeof(MFileChain));
			if (ifs.eof()) { throw std::exception(); }
			chains.emplace_back(chain);

			std::vector< MFileBase1> basetracks;
			basetracks.reserve(chain.nseg);

			for (int b = 0; b < chain.nseg; b++) {
				MFileBase1 base;
				ifs.read((char*)&base, sizeof(MFileBase));
				if (ifs.eof()) { throw std::exception(); }
				basetracks.emplace_back(base);
			}
			all_basetracks.emplace_back(basetracks);
		}
		auto nowpos = ifs.tellg();
		auto size1 = nowpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;

		size_t Nbasetrack = 0;
		for (auto& p : all_basetracks) {
			Nbasetrack += p.size();
		}

		if (Nbasetrack != mfile.info_header.Nbasetrack) { throw std::exception("Nbasetrack is wrong."); }

		mfile.chains.swap(chains);
		mfile.all_basetracks.swap(all_basetracks);
	}
	void mfile1::read_mfile(std::string filepath, MFile_minimum& mfile) {
		std::ifstream ifs(filepath, std::ios::binary);
		//filesizeéÊìæ
		ifs.seekg(0, std::ios::end);
		int64_t eofpos = ifs.tellg();
		ifs.clear();
		ifs.seekg(0, std::ios::beg);
		int64_t begpos = ifs.tellg();
		int64_t size2 = eofpos - begpos;
		int64_t GB = size2 / (1000 * 1000 * 1000);
		int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
		int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
		if (GB > 0) {
			std::cout << "FILE size :" << GB << "." << std::setw(3) << std::setfill('0') << MB << " [GB]" << std::endl;
		}
		else {
			std::cout << "FILE size :" << MB << "." << std::setw(3) << std::setfill('0') << KB << " [MB]" << std::endl;
		}
		//Mfile headerÇÃì«Ç›çûÇ›
		ifs.read((char*)&mfile.header, sizeof(MFileHeader));
		if (ifs.eof()) { throw std::exception(); }
		std::string  filetype = "mfile-a0";
		memcpy((char*)filetype.data(), (char*)&mfile.header.filetype, filetype.size());
		if (filetype != "mfile-a0") { throw std::exception("File format is not mfile-a0."); }

		//mfile info headerÇÃì«Ç›çûÇ›
		ifs.read((char*)&mfile.info_header, sizeof(MFileInfoHeader));
		if (ifs.eof()) { throw std::exception(); }

		if (sizeof(MFileChain) != mfile.info_header.classsize1) { throw std::exception("Classsize1 is wrong."); }
		if (sizeof(MFileBase) != mfile.info_header.classsize2) { throw std::exception("Classsize2 is wrong."); }

		std::vector< MFileChain> chains;
		chains.reserve(mfile.info_header.Nchain);

		std::vector< std::vector< MFileBase>> all_basetracks;
		int count = 0;
		for (uint32_t c = 0; c < mfile.info_header.Nchain; c++) {
			if (count % 100000 == 0) {
				auto nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			count++;
			MFileChain chain;
			ifs.read((char*)&chain, sizeof(MFileChain));
			if (ifs.eof()) { throw std::exception(); }
			chains.emplace_back(chain);

			std::vector< MFileBase> basetracks;
			basetracks.reserve(chain.nseg);

			for (int b = 0; b < chain.nseg; b++) {
				MFileBase base;
				ifs.read((char*)&base, sizeof(MFileBase));
				if (ifs.eof()) { throw std::exception(); }
				basetracks.emplace_back(base);
			}
			all_basetracks.emplace_back(basetracks);
		}
		auto nowpos = ifs.tellg();
		auto size1 = nowpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;

		size_t Nbasetrack = 0;
		for (auto& p : all_basetracks) {
			Nbasetrack += p.size();
		}

		if (Nbasetrack != mfile.info_header.Nbasetrack) { throw std::exception("Nbasetrack is wrong."); }

		mfile.chains.swap(chains);
		mfile.all_basetracks.swap(all_basetracks);

	}

	void mfile1::read_mfile_txt(std::string filepath, MFile& mfile) {

		std::ifstream ifs(filepath);
		//filesizeéÊìæ
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
			std::cout << "FILE size :" << GB << "." << std::setw(3) << std::setfill('0') << MB << " [GB]" << std::endl;
		}
		else {
			std::cout << "FILE size :" << MB << "." << std::setw(3) << std::setfill('0') << KB << " [MB]" << std::endl;
		}

		std::string str;

		mfile0::M_Header head_tmp;
		for (int i = 0; i < 3; i++) {
			std::getline(ifs, str);
			head_tmp.head[i] = str;
		}
		std::getline(ifs, str);
		head_tmp.num_all_plate = stoi(str);

		{
			std::getline(ifs, str);
			auto str_v = StringSplit(str);
			for (int i = 0; i < str_v.size(); i++) {
				head_tmp.all_pos.push_back(stoi(str_v[i]));
			}
		}

		int cnt = 0;
		std::vector< MFileChain1> chains;
		chains.reserve(mfile.info_header.Nchain);

		std::vector< std::vector< MFileBase1>> all_basetracks;

		while (std::getline(ifs, str)) {
			if (cnt % 10000 == 0) {
				nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			cnt++;
			MFileChain1 chain;

			auto strs = StringSplit(str);
			if (strs.size() == 4) {
				chain.chain_id = stoi(strs[0]);
				chain.nseg = stoi(strs[1]);
				chain.pos0 = stoi(strs[2]);
				chain.pos1 = stoi(strs[3]);
				std::vector< MFileBase1> basetracks;

				for (int i = 0; i < chain.nseg; i++) {
					MFileBase1 base;
					//int pos, group_id;
					//uint64_t rawid;
					//int ph;
					//float ax, ay;
					//double x, y, z;

					std::getline(ifs, str);
					strs = StringSplit(str);
					if (strs.size() == 15) {
						base.pos = stoi(strs[0]);
						base.group_id = stoi(strs[1]);
						base.rawid = stoi(strs[2]);
						base.ph = stoi(strs[3]);
						base.ax = stof(strs[4]);
						base.ay = stof(strs[5]);
						base.x = stod(strs[6]);
						base.y = stod(strs[7]);
						base.z = stod(strs[8]);
						base.tmp[0] = 0;
						base.tmp[1] = 0;
						base.tmp[2] = 0;
						base.tmp[3] = 0;
						basetracks.emplace_back(base);
					}
					else {
						fprintf(stderr, "mfile format err! cannot read file\n");
						fprintf(stderr, "base track information 15 -->%d\n", int(strs.size()));
						for (int i = 0; i < strs.size(); i++) {
							printf("%s ", strs[i].c_str());
						}
						printf("\n");
						throw std::exception();
					}
				}
				all_basetracks.emplace_back(basetracks);
				chains.emplace_back(chain);
			}
			else {
				fprintf(stderr, "mfile format err! cannot read file\n");
				fprintf(stderr, "chain information 4 -->%d\n", int(strs.size()));
				for (int i = 0; i < strs.size(); i++) {
					printf("%s ", strs[i].c_str());
				}
				printf("\n");

				throw std::exception();
			}
		}
		auto size1 = eofpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
		if (cnt == 0) {
			fprintf(stderr, "%s no chain!\n", filepath.c_str());
			exit(1);
		}


		mfile.chains.swap(chains);
		mfile.all_basetracks.swap(all_basetracks);
	}

	void mfile1::write_mfile(std::string filepath, MFile& mfile) {
		std::ofstream ofs(filepath, std::ios::binary);

		std::string  filetype = "mfile-a0";
		mfile.header.filetype = 0;
		memcpy((char*)&mfile.header.filetype, filetype.data(), filetype.size());

		mfile.info_header.classsize1 = sizeof(MFileChain);
		mfile.info_header.classsize2 = sizeof(MFileBase);

		mfile.info_header.Nchain = mfile.chains.size();
		assert(mfile.chains.size() == mfile.all_basetracks.size());

		size_t Nbasetrack = 0;
		for (auto& p : mfile.all_basetracks) {
			Nbasetrack += p.size();
		}
		mfile.info_header.Nbasetrack = Nbasetrack;

		ofs.write((char*)&mfile.header, sizeof(MFileHeader));
		ofs.write((char*)&mfile.info_header, sizeof(MFileInfoHeader));

		int count = 0;
		int64_t max = mfile.info_header.Nchain;

		for (int c = 0; c < mfile.info_header.Nchain; c++) {
			if (count % 10000 == 0) {
				std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
			}
			count++;
			ofs.write((char*)&mfile.chains[c], sizeof(MFileChain));
			assert(mfile.chains[c].nseg == mfile.all_basetracks[c].size());
			for (int b = 0; b < mfile.chains[c].nseg; b++) {
				ofs.write((char*)&mfile.all_basetracks[c][b], sizeof(MFileBase));
			}
		}
		std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;

	}
	void mfile1::write_mfile(std::string filepath, MFile_minimum& mfile) {
		std::ofstream ofs(filepath, std::ios::binary);

		std::string  filetype = "mfile-a0";
		mfile.header.filetype = 0;
		memcpy((char*)&mfile.header.filetype, filetype.data(), filetype.size());

		mfile.info_header.classsize1 = sizeof(MFileChain);
		mfile.info_header.classsize2 = sizeof(MFileBase);

		mfile.info_header.Nchain = mfile.chains.size();
		assert(mfile.chains.size() == mfile.all_basetracks.size());

		size_t Nbasetrack = 0;
		for (auto& p : mfile.all_basetracks) {
			Nbasetrack += p.size();
		}
		mfile.info_header.Nbasetrack = Nbasetrack;

		ofs.write((char*)&mfile.header, sizeof(MFileHeader));
		ofs.write((char*)&mfile.info_header, sizeof(MFileInfoHeader));

		int count = 0;
		int64_t max = mfile.info_header.Nchain;

		for (int c = 0; c < mfile.info_header.Nchain; c++) {
			if (count % 10000 == 0) {
				std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
			}
			count++;
			ofs.write((char*)&mfile.chains[c], sizeof(MFileChain));
			assert(mfile.chains[c].nseg == mfile.all_basetracks[c].size());
			for (int b = 0; b < mfile.chains[c].nseg; b++) {
				ofs.write((char*)&mfile.all_basetracks[c][b], sizeof(MFileBase));
			}
		}
		std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;

	}
	void mfile1::read_mfile_extension(std::string filename, mfile0::Mfile& m) {
		m.chains.reserve(100000000);
		std::string extension;
		extension = filename.substr(filename.size() - 3, 3);
		if (extension == "bmf") {
			mfile1::MFile m_tmp;
			mfile1::read_mfile(filename, m_tmp);
			mfile1::converter(m_tmp, m);
		}
		else {
			mfile0::read_mfile(filename, m);
		}
		m.chains.shrink_to_fit();
	}
	void mfile1::write_mfile_extension(std::string filename, mfile0::Mfile& m) {
		std::string extension;
		extension = filename.substr(filename.size() - 3, 3);
		if (extension == "bmf") {
			mfile1::MFile m_tmp;
			mfile1::converter(m, m_tmp);
			mfile1::write_mfile(filename, m_tmp);
		}
		else {
			mfile0::write_mfile(filename, m);
		}
	}



	double mfile1::MFile::chain_ax(int i) {
		double tmp = 0;
		int count = 0;
		for (auto itr = all_basetracks[i].begin(); itr != all_basetracks[i].end(); itr++) {
			tmp += itr->ax;
			count++;
		}
		return tmp / count;
	}
	double mfile1::MFile::chain_ay(int i) {
		double tmp = 0;
		int count = 0;
		for (auto itr = all_basetracks[i].begin(); itr != all_basetracks[i].end(); itr++) {
			tmp += itr->ay;
			count++;
		}
		return tmp / count;
	}
	double mfile1::MFile::angle_diff_dev_rad(int i, double ax, double ay) {
		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double angle = sqrt(ax * ax + ay * ay);
		int count = 0;
		for (auto itr = all_basetracks[i].begin(); itr != all_basetracks[i].end(); itr++) {
			diff = (itr->ax * ax + itr->ay * ay) / angle - angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
		}
		return sqrt(sum2 / count - pow(sum / count, 2.));
	}
	double mfile1::MFile::angle_diff_dev_lat(int i, double ax, double ay) {
		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double angle = sqrt(ax * ax + ay * ay);
		int count = 0;
		for (auto itr = all_basetracks[i].begin(); itr != all_basetracks[i].end(); itr++) {
			diff = (itr->ay * ax - itr->ay * ax) / angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
		}
		return sqrt(sum2 / count - pow(sum / count, 2.));
	}
	double chain_ax(std::vector<MFileBase>& b) {
		double tmp = 0;
		int count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			tmp += itr->ax;
			count++;
		}
		return tmp / count;
	}
	double chain_ay(std::vector<MFileBase>& b) {
		double tmp = 0;
		int count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			tmp += itr->ay;
			count++;
		}
		return tmp / count;

	}
	double angle_diff_dev_rad(std::vector<MFileBase>& b) {
		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double ax = chain_ax(b);
		double ay = chain_ay(b);
		double angle = sqrt(ax * ax + ay * ay);
		int count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			diff = (itr->ax * ax + itr->ay * ay) / angle - angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
		}
		return sqrt(sum2 / count - pow(sum / count, 2.));

	}
	double angle_diff_dev_lat(std::vector<MFileBase>& b) {
		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double ax = chain_ax(b);
		double ay = chain_ay(b);
		double angle = sqrt(ax * ax + ay * ay);
		int count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			diff = (itr->ay * ax - itr->ax * ay) / angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
		}
		return sqrt(sum2 / count - pow(sum / count, 2.));

	}
	double chain_vph(std::vector<MFileBase>& b) {
		double tmp = 0;
		int count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			tmp += itr->ph % 10000;
			count++;
		}
		return tmp / count;
	}
	double chain_ph(std::vector<MFileBase>& b) {
		double tmp = 0;
		int count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			tmp += (int)(itr->ph / 10000);
			count++;
		}
		return tmp / count;
	}
	double angle_diff_mom_iron(std::vector<MFileBase>& b, int num_thr) {

		double sum = 0;
		double sum2 = 0;
		double diff = 0;
		double ax, ay, angle;
		ax = 0;
		ay = 0;
		int pl, count = 0;
		for (auto itr = b.begin(); itr != b.end(); itr++) {
			if (itr + 1 == b.end())continue;
			if ((itr + 1)->pos / 10 - itr->pos / 10 != 1)continue;
			pl = itr->pos / 10;
			if (pl == 1 || pl == 2 || pl == 3 || pl == 15)continue;
			if (pl > 15 && pl % 2 == 1)continue;
			angle = sqrt(itr->ax * itr->ax + itr->ay * itr->ay);
			diff = ((itr + 1)->ay * itr->ax - (itr + 1)->ax * itr->ay) / angle;
			sum += diff;
			sum2 += diff * diff;
			count++;
			ax += ((itr + 1)->ax + itr->ax) / 2;
			ay += ((itr + 1)->ay + itr->ay) / 2;
		}
		if (count < num_thr)return -1;
		ax = ax / count;
		ay = ay / count;
		double rms = sqrt(sum2 / count - pow(sum / count, 2.));
		double iron_thick = 500 * sqrt(1 + ax * ax + ay * ay) * count;
		double iron_rad = 17.57 * 1000;
		double bpc = 13.6 / rms * sqrt(iron_thick / iron_rad) * (1 + 0.038 * std::log(iron_thick / iron_rad));
		return bpc;

	}

}

