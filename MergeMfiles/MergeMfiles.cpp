#pragma comment(lib, "FILE_structure.lib")
#include <FILE_structure.hpp>

#include <iomanip>

class output_format {
public:
	int groupid, chainid, pl, upl;
	int peke, count;
	double dal, dar, md, dz, dz2, oa;
};

struct Basetrack {
	int pos, eid, rawid, ph;
	double ax, ay, x, y;
	int z;
	int flg, b, c, d;
	double p, q;
};

struct kink_cand {
	int pl, groupid;
};
bool operator<(const kink_cand& lhs, const kink_cand& rhs) {
	return std::tie(lhs.groupid, lhs.pl) < std::tie(rhs.groupid, rhs.pl);
}
void SetMfile(std::string input, std::multimap<int, Basetrack>& map);
void Output(std::string output, std::multimap<int, Basetrack>& out);
std::multimap<int, Basetrack>  SetAndMatchMfile(std::string input, std::multimap<int, Basetrack>& map);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage : input.all input2.all output.all\n");
		exit(1);
	}
	std::string in_old = argv[1];// input mfile
	std::string in_new = argv[2];// input mfile
	std::string output = argv[3];// output mfile

	std::multimap<int, Basetrack> trk;
	int mode = 1;
	if (mode == 0) {
		SetMfile(in_old, trk);
		SetMfile(in_new, trk);
	}
	else {
		std::multimap<int, Basetrack> tmp;
		SetMfile(in_old, tmp);
		trk = SetAndMatchMfile(in_new, tmp);
	}

	Output(output, trk);
}

void SetMfile(std::string input, std::multimap<int, Basetrack>& map)
{
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::cout << input << std::endl;
	int count = 0;
	Basetrack b;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);


		if (count < 6) {
			// debag
			if (str_v.size() == 4) {
				std::cout << str << std::endl;
			}
			if (str_v.size() == 15) {
				std::cout << str << std::endl;
			}
		}
		if (str_v.size() == 15) {
			b.pos = std::stoi(str_v[0]);
			b.eid = std::stoi(str_v[1]);
			b.rawid = std::stoi(str_v[2]);
			b.ph = std::stoi(str_v[3]);
			b.ax = std::stod(str_v[4]);
			b.ay = std::stod(str_v[5]);
			b.x = std::stod(str_v[6]);
			b.y = std::stod(str_v[7]);
			b.z = std::stod(str_v[8]);
			b.flg = std::stoi(str_v[9]);
			b.b = std::stoi(str_v[10]);

			b.c = std::stoi(str_v[11]);
			b.d = std::stoi(str_v[12]);
			b.p = std::stod(str_v[13]);
			b.q = std::stod(str_v[14]);

		}

		map.insert(std::make_pair(b.eid,b));
		count++;
	}
	std::cout << "\t# of Starttrack : " << count << std::endl;

}

std::multimap<int, Basetrack>  SetAndMatchMfile(std::string input, std::multimap<int, Basetrack>& map)
{
	std::multimap<int, Basetrack> out;
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::cout << input << std::endl;
	int count = 0;
	Basetrack b;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);


		if (count < 10) {
			// debag
			if (str_v.size() == 4) {
				std::cout << str << std::endl;
			}
			if (str_v.size() == 15) {
				std::cout << str << std::endl;
			}
		}
		if (str_v.size() == 15) {
			if (map.find(std::stoi(str_v[1])) != map.end()) {
				b.pos = std::stoi(str_v[0]);
				b.eid = std::stoi(str_v[1]);
				b.rawid = std::stoi(str_v[2]);
				b.ph = std::stoi(str_v[3]);
				b.ax = std::stod(str_v[4]);
				b.ay = std::stod(str_v[5]);
				b.x = std::stod(str_v[6]);
				b.y = std::stod(str_v[7]);
				b.z = std::stod(str_v[8]);
				b.flg = std::stoi(str_v[9]);
				b.b = std::stoi(str_v[10]);

				b.c = std::stoi(str_v[11]);
				b.d = std::stoi(str_v[12]);
				b.p = std::stod(str_v[13]);
				b.q = std::stod(str_v[14]);
				
				for (auto itr = map.upper_bound(b.eid); itr != map.lower_bound(b.eid); itr++) {
					map.insert(std::make_pair(itr->first, itr->second));
				}
				out.insert(std::make_pair(b.eid, b));

			}

		}

		count++;
	}
	std::cout << "\t# of Starttrack : " << count << std::endl;
	return out;
}

void Output(std::string output, std::multimap<int, Basetrack>& out) 
{

	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "File open error : " << output << std::endl;
		exit(1);
	}
	int count = 0;
	
	//header
	ofs << "% Created by mkmf" << std::endl;
	ofs << "% on 2007 / 12 / 25 23:28 : 09 + 09:00 (JST)" << std::endl;
	ofs << "0       0   3   0      0.0   0.0000" << std::endl;
	ofs << "133" << std::endl;
	ofs << "11   21   31   41   51   61   71   81   91  101  111  121  131  141  151  161  171  181  191  201  211  221  231  241  251  261  271  281  291  301  311  321  331  341  351  361  371  381  391  401  411  421  431  441  451  461  471  481  491  501  511  521  531  541  551  561  571  581  591  601  611  621  631  641  651  661  671  681  691  701  711  721  731  741  751  761  771  781  791  801  811  821  831  841  851  861  871  881  891  901  911  921  931  941  951  961  971  981  991 1001 1011 1021 1031 1041 1051 1061 1071 1081 1091 1101 1111 1121 1131 1141 1151 1161 1171 1181 1191 1201 1211 1221 1231 1241 1251 1261 1271 1281 1291 1301 1311 1321 1331" << std::endl;

	std::set<int>set;
	for (auto itr = out.begin(); itr != out.end(); itr++) {
		set.insert(itr->first);
	}

	for (auto itr = set.begin(); itr != set.end(); itr++) {
		
		auto p = out.equal_range(*itr);
		if (p.second != out.end()) {
			ofs << std::fixed << std::right
				<< std::setw(10) << std::setprecision(0) << 0 << " "
				<< std::setw(5) << std::setprecision(0) << out.count(*itr) << " "
				<< std::setw(6) << std::setprecision(0) << p.first->second.pos << " "
				<< std::setw(6) << std::setprecision(0) << std::prev(p.second)->second.pos
				<< std::endl;

			for (auto q = p.first; q != p.second; q++) {
				ofs << std::fixed << std::right
					<< std::setw(6) << std::setprecision(0) << q->second.pos << " "
					<< std::setw(6) << std::setprecision(0) << q->second.eid << " "
					<< std::setw(10) << std::setprecision(0) << q->second.rawid << " "
					<< std::setw(8) << std::setprecision(0) << q->second.ph << " "
					<< std::setw(10) << std::setprecision(4) << q->second.ax << " "
					<< std::setw(10) << std::setprecision(4) << q->second.ay << " "
					<< std::setw(10) << std::setprecision(1) << q->second.x << " "
					<< std::setw(10) << std::setprecision(1) << q->second.y << " "
					<< std::setw(8) << std::setprecision(0) << q->second.z << " "
					<< std::setw(3) << std::setprecision(0) << q->second.flg << " "
					<< std::setw(3) << std::setprecision(0) << q->second.b << " "
					<< std::setw(3) << std::setprecision(0) << q->second.c << " "
					<< std::setw(3) << std::setprecision(0) << q->second.d << " "
					<< std::setw(8) << std::setprecision(4) << q->second.p << " "
					<< std::setw(8) << std::setprecision(4) << q->second.q
					<< std::endl;
			}
		}
		count++;
	}

	std::cout << "writing file is finished!" << std::endl;
}