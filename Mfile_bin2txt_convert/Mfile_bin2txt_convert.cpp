//20231209
//to convert mfile.bin to mfile txt
//kasumi

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.h>
#include <functions.h>
#include <set>
#include <random>

mfile1::MFile mfile_read_rand(std::string file_in_mfile, std::set<uint64_t>& read_gid);
uint64_t get_rand_range(uint64_t min_val, uint64_t max_val);
void read_chain(std::ifstream& ifs, mfile1::MFileChain1& chain, std::vector< mfile1::MFileBase1>& basetracks);

int main(int argc, char** argv) {
	if (argc != 3) {
		//fprintf(stderr, "usage : prg_name [input m-file-bin] groupnum gmin gmax [output m-file-bin]\n");
		fprintf(stderr, "usage : prg_name [input m-file-bin] [output m-file-txt]\n");
		exit(1);
	}
	std::string file_in_mfile = argv[1];
	//uint64_t gnum = std::stoll(argv[2]);
	//uint64_t gmin = std::stoll(argv[3]);
	//uint64_t gmax = std::stoll(argv[4]);
	std::string file_out_mfile = argv[2];

	//std::set<uint64_t>gid;
	//for (uint64_t i = 0; i < gnum; i++) {
	//	gid.insert(get_rand_range(gmin, gmax));
	//}
	//mfile1::MFile m = mfile_read_rand(file_in_mfile, gid);
	
	//binaryの読み込みとmへの突っ込み 
	mfile1::MFile m;
	mfile1::read_mfile(file_in_mfile, m);



	mfile0::Mfile m_txt;
	mfile1::converter(m, m_txt);

	mfile0::write_mfile(file_out_mfile, m_txt);

}
mfile1::MFile mfile_read_rand(std::string file_in_mfile, std::set<uint64_t>& read_gid) {
	std::ifstream ifs(file_in_mfile, std::ios::binary);

	//filesize取得
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
	mfile1::MFile mfile;
	//Mfile headerの読み込み
	ifs.read((char*)&mfile.header, sizeof(mfile1::MFileHeader));
	if (ifs.eof()) { throw std::exception(); }
	std::string  filetype = "mfile-a0";
	memcpy((char*)filetype.data(), (char*)&mfile.header.filetype, filetype.size());
	if (filetype != "mfile-a0") { throw std::exception("File format is not mfile-a0."); }

	//mfile info headerの読み込み
	ifs.read((char*)&mfile.info_header, sizeof(mfile1::MFileInfoHeader));
	if (ifs.eof()) { throw std::exception(); }

	if (sizeof(mfile1::MFileChain) != mfile.info_header.classsize1) { throw std::exception("Classsize1 is wrong."); }
	if (sizeof(mfile1::MFileBase) != mfile.info_header.classsize2) { throw std::exception("Classsize2 is wrong."); }
	//Mfile headerの書き込み
	//mfile info headerの書き込み

	std::vector< mfile1::MFileChain1> group;
	std::vector< std::vector< mfile1::MFileBase1>> all_basetracks;
	uint64_t count = 0, r_base_num = 0, r_chain_num = 0, r_group_num = 0;
	int64_t gid = -1;
	int num;
	for (int64_t c = 0; c < mfile.info_header.Nchain; c++) {
		if (count % 100000 == 0) {
			auto nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;
		mfile1::MFileChain1 chain;
		std::vector< mfile1::MFileBase1> basetracks;
		//chainの読み込み
		read_chain(ifs, chain, basetracks);

		//読んだchainが前のchainと同じgroupか確認
		//最後のchainを読み込んだ時
		if (c + 1 == mfile.info_header.Nchain) {
			//最後のchainが違うgroups
			if (gid != basetracks.begin()->group_id) {
				////////////////////////////
				//今までのgroupの書き出し
				if (read_gid.count(r_group_num) == 1) {
					for (auto c : group) {
						mfile.chains.push_back(c);
					}
					for (auto b : all_basetracks) {
						mfile.all_basetracks.push_back(b);
					}
				}
				r_group_num++;
				gid = basetracks.begin()->group_id;
				group.clear();
				for (int64_t i = 0; i < all_basetracks.size(); i++) {
					all_basetracks[i].clear();
				}
				all_basetracks.clear();
				////////////////////////////
				//最後のgroupの書き出し
				all_basetracks.emplace_back(basetracks);
				group.emplace_back(chain);
				if (read_gid.count(r_group_num) == 1) {
					for (auto c : group) {
						mfile.chains.push_back(c);
					}
					for (auto b : all_basetracks) {
						mfile.all_basetracks.push_back(b);
					}
				}
				r_group_num++;
				gid = basetracks.begin()->group_id;
				group.clear();
				for (int64_t i = 0; i < all_basetracks.size(); i++) {
					all_basetracks[i].clear();
				}
				all_basetracks.clear();
			}
			else {
				//今のchainをpush back
				all_basetracks.emplace_back(basetracks);
				group.emplace_back(chain);
				if (read_gid.count(r_group_num) == 1) {
					for (auto c : group) {
						mfile.chains.push_back(c);
					}
					for (auto b : all_basetracks) {
						mfile.all_basetracks.push_back(b);
					}
				}
				r_group_num++;
				gid = basetracks.begin()->group_id;
				////////////////////////////
				//今までのgroupの書き出し
				r_group_num++;
				gid = basetracks.begin()->group_id;
				group.clear();
				for (int64_t i = 0; i < all_basetracks.size(); i++) {
					all_basetracks[i].clear();
				}
				all_basetracks.clear();
			}
		}
		else if (group.size() != 0 && gid != basetracks.begin()->group_id) {
			if (read_gid.count(r_group_num) == 1) {
				for (auto c : group) {
					mfile.chains.push_back(c);
				}
				for (auto b : all_basetracks) {
					mfile.all_basetracks.push_back(b);
				}
			}
			r_group_num++;
			gid = basetracks.begin()->group_id;
			group.clear();
			for (int64_t i = 0; i < all_basetracks.size(); i++) {
				all_basetracks[i].clear();
			}
			all_basetracks.clear();
		}
		else if (c == 0) {
			gid = basetracks.begin()->group_id;
		}
		all_basetracks.emplace_back(basetracks);
		group.emplace_back(chain);

	}
	auto nowpos = ifs.tellg();
	auto size1 = nowpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;

	return mfile;
}
void read_chain(std::ifstream& ifs, mfile1::MFileChain1& chain, std::vector< mfile1::MFileBase1>& basetracks) {

	//一旦chainを読む
	ifs.read((char*)&chain, sizeof(mfile1::MFileChain));
	if (ifs.eof()) { throw std::exception(); }
	basetracks.clear();
	basetracks.reserve(chain.nseg);

	for (int b = 0; b < chain.nseg; b++) {
		mfile1::MFileBase1 base;
		ifs.read((char*)&base, sizeof(mfile1::MFileBase));
		if (ifs.eof()) { throw std::exception(); }
		basetracks.emplace_back(base);
	}
}
uint64_t get_rand_range(uint64_t min_val, uint64_t max_val) {
	// 乱数生成器
	static std::mt19937_64 mt64(0);

	// [min_val, max_val] の一様分布整数 (int) の分布生成器
	std::uniform_int_distribution<uint64_t> get_rand_uni_int(min_val, max_val);

	// 乱数を生成
	return get_rand_uni_int(mt64);
}