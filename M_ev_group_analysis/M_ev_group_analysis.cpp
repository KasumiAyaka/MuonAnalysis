//calc_some_value_of_indicated_events
// 2025/01/22
// 2024/09/07
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>

class output_format {
public:
	int groupid, chainid, pl, upl;
	int peke, count;
	double dal, dar, md, dz, dz2, oa;
};


struct Key {
	int pl, cid, gid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.pl) < std::tie(rhs.gid, rhs.cid, rhs.pl);
}

struct Basetrack {
	Key ky;
	int vph, rid;
	double x, y, z, ax, ay;
};
struct Btrk {
	Basetrack b0, b1;
};
void Partner_analysis(mfile0::Mfile& m_all, mfile0::Mfile& m_all_out);


int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:file-in.all output.all\n");
		exit(1);
	}
	std::string file_in_mfile = argv[1];
	std::string file_out_mfile = argv[2];
	//std::string file_in_list = argv[2];
	//int eccnum = std::stoi(argv[3]);
	//std::string file_out = argv[4];

	mfile0::Mfile all, out_all;
	mfile1::read_mfile_extension(file_in_mfile, all);

	Partner_analysis(all, out_all);

	mfile1::write_mfile_extension(file_out_mfile, out_all);
}

void Partner_analysis(mfile0::Mfile& m_all, mfile0::Mfile& m_all_out) {
	std::pair<int64_t, int> t,t0,t1;
	std::multimap<std::tuple<int64_t, int>, mfile0::M_Chain> map;
	std::set<int64_t> set;

	Btrk tmp;

	int i = 0;
	int64_t gid = 0;
	int flg = 0;

	for (auto itr = m_all.chains.begin(); itr != m_all.chains.end(); itr++) {
		gid = itr->basetracks.begin()->group_id;
		//std::cout << "\ngid = " << gid << std::endl;
		flg = gid % 2;
		//std::cout << "flg = " << flg << std::endl;
		gid = gid / 10;
		//std::cout << "gid = " << gid << ", cid = " << itr->chain_id << std::endl;

		t= std::make_pair(gid,flg);
		map.insert(std::make_pair(t, *itr));
		set.insert(gid);
		i++;
		//if (i > 1000)break;
	}
	std::cout << std::endl;



	std::map<int64_t, int> muonpl;
	m_all_out.header = m_all.header;
	int e0, e1;
	int a0, a1;
	for (auto itr = set.begin(); itr != set.end(); itr++) {
		//gid,flg
		t0 = std::make_pair(*itr, 0);
		t1 = std::make_pair(*itr, 1);
		// # chain
		e0 = map.count(t0);
		e1 = map.count(t1);

		//std::cout << "  # of Track (gid, cid, event, partner-check) = ("
		//	<< map.upper_bound(t0)->second.basetracks.begin()->group_id << ", " << map.upper_bound(t0)->second.chain_id << ", " 
		//	<< e0 << ", " << e1 << ")" << std::endl;

		auto p = map.equal_range(t0);
		auto q = map.equal_range(t1);
		int flg = 0;
		double mu_x, mu_y;
		if (e1 - e0 > 20 || e1 - e0 < 1) {
			// partner‚Ìchain convolusion‚ª‰˜‚¢‚â‚Â‚Æpartner‚ÉŒó•â‚ª‚¢‚È‚¢‚â‚Â
			//std::cout << "  # of Track (gid, cid, event, partner-check) = ("
			//	<< p.first->second.basetracks.begin()->group_id << ", " << p.first->second.chain_id << ", "
			//	<< e0 << ", " << e1 << ")" << std::endl;
		}
		/*
		else {
			for (auto itr0 = p.first; itr0 != p.second; itr0++) {
				m_all_out.chains.push_back(itr0->second);
			}

			for (auto itr0 = q.first; itr0 != q.second; itr0++) {
				m_all_out.chains.push_back(itr0->second);
			}
		}
		*/
		else {
			for (auto itr0 = p.first; itr0 != p.second; itr0++) {
				if (itr0->second.chain_id == 0) {
					mu_x = itr0->second.basetracks.rbegin()->x;
					mu_y = itr0->second.basetracks.rbegin()->y;
					if (mu_x > 40000 && mu_x < 210000 && mu_y>40000 && mu_y < 210000) {
						flg++;
						std::cout << "cid = " << itr0->second.chain_id << " " << itr0->second.basetracks.rbegin()->pos << std::endl;
					}
				}
				if (flg > 0) {
					m_all_out.chains.push_back(itr0->second);
				}
			}

			for (auto itr0 = q.first; itr0 != q.second; itr0++) {
				if (flg > 0) {
					m_all_out.chains.push_back(itr0->second);
				}
			}
		}

		/*
		else {
			int mu_pl = 0;
			flg = 0;
			for (auto itr0 = p.first; itr0 != p.second; itr0++) {
				std::cout << " tnum==2?" << e0 << std::endl;

				m_all_out.chains.push_back(itr0->second);

				if (itr0->second.chain_id == 0) {
					mu_pl = itr0->second.pos1/10;
					muonpl.insert(std::make_pair(t0.first, mu_pl));
				}
				if (itr0->second.chain_id == 1) {
					a1 = itr0->second.pos1 / 10;
					a0 = itr0->second.pos0 / 10;
					if (a1 > mu_pl && a0 < mu_pl + 1) {
						//partner‚Ì“_‚ÅŠÑ’Ê
						flg = 1;
					}
				}
			}

			for (auto itr0 = q.first; itr0 != q.second; itr0++) {
				if (itr0->second.chain_id <2) {
					m_all_out.chains.push_back(itr0->second);
					std::cout << "\tcheck " << itr0->second.chain_id << std::endl;
				}
				else {
					if (flg > 0) {
						std::encl;
					}
					
				}
			}
		}
		*/
	}
	/*
	partner‚ªmuon stop pl‚ğ‚Ü‚½‚®‚±‚Æ‚ğ—v‹‚·‚é‘O‚Ìrequest:partner‚Ì‰˜‚¢‚â‚Â‚ğœ‹‚·‚é‚¾‚¯
	for (auto itr = set.begin(); itr != set.end(); itr++) {
		t0 = std::make_pair(*itr, 0);
		t1 = std::make_pair(*itr, 1);
		e0 = map.count(t0);
		e1 = map.count(t1);

		//std::cout << "  # of Track (gid, cid, event, partner-check) = ("
		//	<< map.upper_bound(t0)->second.basetracks.begin()->group_id << ", " << map.upper_bound(t0)->second.chain_id << ", " 
		//	<< e0 << ", " << e1 << ")" << std::endl;

		auto p = map.equal_range(t0);
		auto q = map.equal_range(t1);
		if (e1 - e0 > 20) {
			// partner‚Ìchain convolusion‚ª‰˜‚¢‚â‚Â
			std::cout << "  # of Track (gid, cid, event, partner-check) = ("
				<< p.first->second.basetracks.begin()->group_id << ", " << p.first->second.chain_id << ", "
				<< e0 << ", " << e1 << ")" << std::endl;
		}
		else {
			for (auto itr0 = p.first; itr0 != p.second; itr0++) {
				m_all_out.chains.push_back(itr0->second);
			}

			for (auto itr0 = q.first; itr0 != q.second; itr0++) {
				m_all_out.chains.push_back(itr0->second);
			}
		}
	}
	*/

}

