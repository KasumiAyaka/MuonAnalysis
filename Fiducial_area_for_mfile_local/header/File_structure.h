#pragma once

#include <fstream>
#include <iostream>
#include <ios>
#include <iomanip> 
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <map>

#include <list>
#include <set>
#include <unordered_map>
#include <utility>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <cassert>

#include "VxxReader.h"
#include <functions.h>

namespace {
	std::vector<std::string> StringSplit(std::string str) {
		std::stringstream ss{ str };
		std::vector<std::string> v;
		std::string buf;
		while (std::getline(ss, buf, ' ')) {
			if (buf != "") {
				v.push_back(buf);
			}
		}
		return v;
	}
	std::vector<std::string> StringSplit_with_tab(std::string str) {
		std::vector<std::string> v;

		std::vector<std::string> str_v = StringSplit(str);
		for (int i = 0; i < str_v.size(); i++) {
			std::stringstream ss{ str_v[i] };
			std::string buf;
			while (std::getline(ss, buf, '\t')) {
				if (buf != "") {
					v.push_back(buf);
				}
			}
		}
		return v;
	}
}
namespace delaunay {

	struct Point {
		/*** Point on the 2D-plane ***/
		double x, y;
		Point operator + (const Point& p) const { return Point{ x + p.x, y + p.y }; }
		Point operator - (const Point& p) const { return Point{ x - p.x, y - p.y }; }
		Point operator * (double d) const { return Point{ x * d, y * d }; }
		Point operator / (double d) const { return Point{ x / d, y / d }; }
		bool operator < (const Point& p) const { return x == p.x ? y < p.y : x < p.x; }
		Point normal() const { return Point{ y, -x }; }
		double norm() const { return x * x + y * y; }
	};

	double dot(Point a, Point b);
	double cross(Point a, Point b);

	/*** Edge connecting two points ***/
	using Edge = std::pair<size_t, size_t>;
	Edge make_edge(size_t a, size_t b);

	/*** Core Implementation ***/
	class DelaunayTriangulation {
		bool is_ccw(size_t a, size_t b, size_t c) const { return cross(P[b] - P[a], P[c] - P[a]) > 0; }

		struct Triangle {
			/*** Triangle on the 2D-plane ***/
			size_t a, b, c;
			size_t opposite(Edge e) const {
				if (e.first != a and e.second != a) return a;
				if (e.first != b and e.second != b) return b;
				return c;
			}
		};

		Triangle make_triangle(size_t a, size_t b, size_t c) const {
			/*** Make triangle where the direction 'a->b->c' is counter-clockwise ***/
			if (!is_ccw(a, b, c)) std::swap(b, c);
			return Triangle{ a, b, c };
		}

		bool point_in_triangle(size_t tar, const Triangle& t) const {
			/*** Check wheter a point lies in an triangle or not ***/
			if (!is_ccw(t.a, t.b, tar)) return false;
			if (!is_ccw(t.b, t.c, tar)) return false;
			if (!is_ccw(t.c, t.a, tar)) return false;
			return true;
		}

	private:
		size_t n;
		std::vector<Point> P;
		std::vector<Triangle> T;
		std::vector<std::list<size_t>> CH;
		std::vector<Edge> edge;

		size_t add_child_triangle(size_t k, const Triangle& t) {
			/*** Add a child triangle for T[k] ***/
			size_t new_k = T.size();
			T.push_back(t);
			CH.push_back(std::list<size_t>());
			CH[k].push_back(new_k);
			return new_k;
		}

		size_t add_child_triangle(size_t k1, size_t k2, const Triangle& t) {
			/*** Add a child triangle for T[k1] and T[k2] ***/
			size_t new_k = T.size();
			T.push_back(t);
			CH.push_back(std::list<size_t>());
			CH[k1].push_back(new_k);
			CH[k2].push_back(new_k);
			return new_k;
		}

		size_t find_child(size_t k, size_t tar) const {
			/*** Find child triangle of T[k] where P[tar] lies in ***/
			for (size_t i : CH[k]) if (point_in_triangle(tar, T[i])) return i;
			return std::numeric_limits<size_t>::max();
		}

	private:
		struct EdgeHash {
			/*** Hash function for edges ***/
			size_t operator () (const Edge& e) const {
				static const size_t FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
				size_t key = ((e.first << 16) ^ e.second) ^ FIXED_RANDOM;
				return key ^ (key >> 16);
			}
		};

		// Maps an edge to its incident triangles
		std::unordered_map<Edge, std::set<size_t>, EdgeHash> e2t;

		void register_triangle(size_t k, const Triangle& t) {
			/*** Register a new triangle ***/
			e2t[make_edge(t.a, t.b)].insert(k);
			e2t[make_edge(t.b, t.c)].insert(k);
			e2t[make_edge(t.c, t.a)].insert(k);
		}

		void unregister_triangle(size_t k, const Triangle& t) {
			/*** Unregister an old triangle ***/
			e2t[make_edge(t.a, t.b)].erase(k);
			e2t[make_edge(t.b, t.c)].erase(k);
			e2t[make_edge(t.c, t.a)].erase(k);
		}

	private:
		bool is_ilegal(Edge e, size_t c, size_t tar) {
			/*** Check wheter an edge is ilegal or not ***/
			size_t a = e.first, b = e.second;

			Point s1 = (P[a] + P[b]) / 2.;
			Point s2 = (P[b] + P[c]) / 2.;
			Point t1 = s1 + (P[b] - P[a]).normal();
			Point t2 = s2 + (P[c] - P[b]).normal();

			double d1 = cross(t2 - s2, s2 - s1);
			double d2 = cross(t2 - s2, t1 - s1);
			Point ct = s1 + (t1 - s1) * d1 / d2;

			double dx1 = P[a].x - ct.x;
			double dy1 = P[a].y - ct.y;
			double dx2 = P[tar].x - ct.x;
			double dy2 = P[tar].y - ct.y;
			return (dx1 * dx1 + dy1 * dy1) > (dx2 * dx2 + dy2 * dy2);
		}

		void legalize_edge(size_t piv, Edge e) {
			/*** Legalize an edge with pivot P[piv] ***/
			if (e2t.count(e) == 0) return;
			if (e2t[e].size() != 2) return;

			size_t k1 = *e2t[e].begin();
			size_t k2 = *e2t[e].rbegin();
			size_t a = T[k1].opposite(e);
			size_t b = T[k2].opposite(e);

			if (is_ilegal(e, a, b)) {
				unregister_triangle(k1, T[k1]);
				unregister_triangle(k2, T[k2]);
				e2t.erase(e);

				Triangle t1 = make_triangle(e.first, a, b);
				Triangle t2 = make_triangle(e.second, a, b);

				size_t k1_ = add_child_triangle(k1, k2, t1);
				size_t k2_ = add_child_triangle(k1, k2, t2);

				register_triangle(k1_, t1);
				register_triangle(k2_, t2);

				size_t c = (piv != a ? a : b);
				legalize_edge(piv, make_edge(e.first, c));
				legalize_edge(piv, make_edge(e.second, c));
			}
		}

		void sub_division(size_t k, size_t tar) {
			/*** Divide triangle T[k] by P[tar] ***/
			unregister_triangle(k, T[k]);

			Triangle t1 = make_triangle(T[k].a, T[k].b, tar);
			Triangle t2 = make_triangle(T[k].b, T[k].c, tar);
			Triangle t3 = make_triangle(T[k].c, T[k].a, tar);

			size_t k1 = add_child_triangle(k, t1);
			size_t k2 = add_child_triangle(k, t2);
			size_t k3 = add_child_triangle(k, t3);

			register_triangle(k1, t1);
			register_triangle(k2, t2);
			register_triangle(k3, t3);

			legalize_edge(tar, make_edge(T[k].a, T[k].b));
			legalize_edge(tar, make_edge(T[k].b, T[k].c));
			legalize_edge(tar, make_edge(T[k].c, T[k].a));
		}

		void modify_convexity() {
			/*** Modify edges of the bounding convex hull ***/
			for (size_t a = n; a <= n + 2; ++a) {
				for (size_t b = 0; b < n; ++b) {
					Edge e(b, a);
					if (e2t.count(e) == 0) continue;

					size_t k1 = *e2t[e].begin();
					size_t k2 = *e2t[e].rbegin();

					size_t c = T[k1].opposite(e);
					size_t d = T[k2].opposite(e);

					if (!is_ccw(a, b, c)) std::swap(c, d);
					if (!is_ccw(c, d, b)) continue;

					unregister_triangle(k1, T[k1]);
					unregister_triangle(k2, T[k2]);

					Triangle t1 = make_triangle(a, c, d);
					Triangle t2 = make_triangle(b, c, d);

					size_t k1_ = add_child_triangle(k1, k2, t1);
					size_t k2_ = add_child_triangle(k1, k2, t2);

					register_triangle(k1_, t1);
					register_triangle(k2_, t2);
				}
			}
		}

	private:
		std::uint32_t seed;

		std::uint32_t xor128() {
			static std::uint32_t x = seed, y = 362436069, z = 521288629, w = 88675123;
			std::uint32_t t = (x ^ (x << 11));
			x = y; y = z; z = w;
			return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
		}

	public:
		DelaunayTriangulation(const std::vector<double>& x, const std::vector<double>& y, uint32_t seed_ = 123456789) {
			n = x.size();
			P = std::vector<Point>(n);
			for (size_t i = 0; i < n; ++i) {
				P[i].x = x[i];
				P[i].y = y[i];
			}

			double r = std::numeric_limits<double>::min();
			for (size_t i = 0; i < n; ++i) {
				r = std::max(r, std::abs(P[i].x));
				r = std::max(r, std::abs(P[i].y));
			}

			// Add three vertices of the initial bounding triangle
			P.push_back(Point{ 3.1 * r, 0 });
			P.push_back(Point{ 0, 3.1 * r });
			P.push_back(Point{ -3.1 * r, -3.1 * r });

			seed = seed_;
		}

		void execute(double min_delta = 1e-6, double max_delta = 1e-5, int max_miss_count = 30) {
			/*** Execute the algorithm ***/

			// Random reordering
			std::vector<size_t> id(n);
			std::iota(id.begin(), id.end(), size_t(0));
			for (size_t i = 0; i < n; ++i) {
				size_t j = xor128() % (n - i);
				std::swap(id[j], id[n - i - 1]);
			}

			// Initialize each data structure
			T = std::vector<Triangle>();
			CH = std::vector<std::list<size_t>>();
			T.push_back(make_triangle(n, n + 1, n + 2));
			CH.push_back(std::list<size_t>());
			register_triangle(0, T[0]);
			e2t = std::unordered_map<Edge, std::set<size_t>, EdgeHash>();

			// Main iteration
			for (size_t tar : id) {
				int miss_count = 0;
				double px = P[tar].x, py = P[tar].y;
				size_t k = 0;
				while (!CH[k].empty()) {
					size_t nxt = find_child(k, tar);
					if (nxt == std::numeric_limits<size_t>::max()) {
						++miss_count;
						if (miss_count >= max_miss_count) break;
						// P[tar] perturbates when it falls on an edge
						double dx = min_delta + (max_delta - min_delta) * xor128() / 0xFFFFFFFFu;
						double dy = min_delta + (max_delta - min_delta) * xor128() / 0xFFFFFFFFu;
						dx *= (xor128() % 2 == 0 ? 1 : -1);
						dy *= (xor128() % 2 == 0 ? 1 : -1);
						P[tar].x = px + dx;
						P[tar].y = py + dy;
					}
					else {
						k = nxt;
					}
				}
				if (CH[k].empty()) sub_division(k, tar);
				P[tar].x = px;
				P[tar].y = py;
			}

			// Modify edges of the bounding convex hull
			modify_convexity();

			// Save the solution
			edge = std::vector<Edge>();
			for (auto it : e2t) {
				Edge e = it.first;
				if (e.first < n and e.second < n) edge.push_back(e);
			}
		}

		const std::vector<Edge>& get_edges() { return edge; }

		void dump(std::ostream& os) {
			os << edge.size() << std::endl;
			for (Edge e : edge) os << e.first << " " << e.second << std::endl;
		}
	};

}

void GaussJorden(double in[3][3], double b[3], double c[3]);

namespace corrmap0 {
	struct Corrmap {
		int id, pos[2];
		double areax[2], areay[2], position[6], angle[6], dz, signal, background, SN, rms_pos[2], rms_angle[2];
		double notuse_d[5];
		int notuse_i[9];
	};
	struct CorrmapDC {
		int id, pos[2], ix, iy;
		double areax[2], areay[2], shr, ddz, dax, day, dz, signal, background, SN, rms_pos[2], rms_angle[2];
		double notuse_d[12];
		int notuse_i[7];
	};
	bool sort_corrid(const Corrmap& left, const Corrmap& right);
	bool sort_corrx(const Corrmap& left, const Corrmap& right);
	void read_cormap(std::string file, std::vector<struct Corrmap>& cor, int output = 1);
	void read_cormap(std::string file, Corrmap& cor);
	void write_corrmap(std::string file, std::vector<struct Corrmap> cor);
	void write_corrmap(std::string file, struct Corrmap cor);
	void read_cormap_dc(std::string file, std::vector<struct CorrmapDC>& cor, int output = 1);

}

//Mfile txt
namespace mfile0 {
	struct M_Header {
		std::string head[3];
		int num_all_plate;
		std::vector<int> all_pos;
	};
	struct M_Base {
		int64_t  rawid, group_id;
		int pos, ph, flg_i[4];
		double x, y, z, flg_d[2];
		float ax, ay;
	};
	struct M_Chain {
		int64_t chain_id;
		int  nseg, pos0, pos1;
		std::vector<M_Base> basetracks;
	};
	struct Mfile {
		M_Header header;
		std::vector<M_Chain> chains;
	};
	void read_mfile(std::string file, Mfile& mfile);
	void read_mfile(std::string file_path, Mfile& mfile, int nseg_thr);
	void write_mfile(std::string filename, const Mfile& mfile, int output = 1);
	void write_mfile_header(std::ofstream& ofs, const M_Header& header, int output = 1);
	void write_mfile_chain(std::ofstream& ofs, const M_Chain& chains);
	double chain_ax(M_Chain chain);
	double chain_ay(M_Chain chain);
	double angle_diff_dev_rad(M_Chain chain, double ax, double ay);
	double angle_diff_dev_lat(M_Chain chain, double ax, double ay);
	double chain_vph(mfile0::M_Chain chain);
	void angle_diff_dev_theta(M_Chain& chain, double& ax, double& ay, double& div_rad, double& div_lat);

	void apply_vph_correction(mfile0::Mfile& m, std::string filename_corr);
	void write_mfile_chain_IVE(std::ofstream& ofs, const M_Chain& chains);
	double chain_fit_dist(M_Chain chain, double slope[3], double intercept[3], double& dist);

	struct M_Chain_inf {
		int64_t chainID;
		int  nseg, pos0, pos1, groupID;
		double ax, ay, radial_deviation, lateral_deviation, d_radial_deviation, d_lateral_deviation, vph_ratio, vph_slope, vph_slope_acc;
		std::pair<double, double> x_up, y_up, x_down, y_down;
		bool outflg_up, outflg_down;
	};
	M_Chain_inf chain2inf(mfile0::M_Chain c);
	void chain_inf_flg(M_Chain_inf& inf, mfile0::M_Chain c, int pl0, int pl1, double range[4], std::map<int, double> z);
	void write_chaininf(std::string filename, std::vector<M_Chain_inf> chain_inf);

	void set_header(int pl0, int pl1, mfile0::Mfile& m);
}
//Mfile bin
namespace mfile1 {
	struct MFileHeader
	{
		uint64_t filetype;
		uint64_t filesize;
		uint64_t reserved;
		uint32_t offset1;
		uint32_t offset2;
		uint8_t name[64];
		uint8_t object[256];
	};

	struct MFileInfoHeader
	{
		uint16_t classsize1;
		uint16_t classsize2;
		uint32_t reserved1;
		uint64_t Nchain;
		uint64_t Nbasetrack;
		uint64_t reserved2;
	};

	class MFileBase {
	public:
		int pos, group_id;
		uint64_t rawid;
		int ph;
		float ax, ay;
		double x, y, z;
	};
	class MFileBase1 : public MFileBase
	{
	public:
		int16_t tmp[4]; // 計算で使うが、ファイルに書く必要のない一時パラメータ
	};


	class MFileChain {
	public:
		uint64_t chain_id;
		int nseg, pos0, pos1;
	};

	class MFileChain1 : public MFileChain
	{
	public:
		int16_t tmp[4]; // 計算で使うが、ファイルに書く必要のない一時パラメータ
	};

	class MFile {

	public:
		MFileHeader header;
		MFileInfoHeader info_header;
		std::vector<MFileChain1> chains;
		std::vector<std::vector<MFileBase1>> all_basetracks;
	public:
		//i番目のchain
		double chain_ax(int i);
		double chain_ay(int i);
		double angle_diff_dev_rad(int i, double ax, double ay);
		double angle_diff_dev_lat(int i, double ax, double ay);


	};
	class MFile_minimum {

	public:
		MFileHeader header;
		MFileInfoHeader info_header;
		std::vector<MFileChain> chains;
		std::vector<std::vector<MFileBase>> all_basetracks;
	public:
		//i番目のchain
		double chain_ax(int i);
		double chain_ay(int i);
		double angle_diff_dev_rad(int i, double ax, double ay);
		double angle_diff_dev_lat(int i, double ax, double ay);

	};
	void converter(const MFile& old, mfile0::Mfile& mfile);
	void converter(const MFile_minimum& old, mfile0::Mfile& mfile);
	void converter(const mfile0::Mfile& old, MFile& mfile);
	void read_mfile(std::string filepath, MFile& mfile);
	void read_mfile(std::string filepath, MFile_minimum& mfile);
	void read_mfile(std::string filepath, MFile_minimum& mfile, int nseg_threshold);
	void read_mfile_txt(std::string filepath, MFile& mfile);

	void write_mfile(std::string filepath, MFile& mfile);
	void write_mfile(std::string filepath, MFile_minimum& mfile);

	void read_mfile_extension(std::string filename, mfile0::Mfile& m);
	void write_mfile_extension(std::string filename, mfile0::Mfile& m);

	double chain_ax(std::vector<MFileBase>& b);
	double chain_ay(std::vector<MFileBase>& b);
	double angle_diff_dev_rad(std::vector<MFileBase>& b);
	double angle_diff_dev_lat(std::vector<MFileBase>& b);
	void angle_diff_dev_theta(std::vector<MFileBase>& b, double& ax, double& ay, double& div_rad, double& div_lat);
	double chain_vph(std::vector<MFileBase>& b);
	double chain_ph(std::vector<MFileBase>& b);
	double angle_diff_mom_iron(std::vector<MFileBase>& b, int num_thr = 5);

}
//vtxfile
namespace vtx {
	class VTX_base {
	public:
		int64_t chainid;
		int pl, rawid, ph;
		double ax, ay, x, y, startPL, nseg, ip;
	};
	class VTX_mindis {
	public:
		int rawid[2], flg;
		double openang, mindis, x, y, z;
	};
	class VTX {
	public:
		int eventid, pl, nbase;
		double x, y, z;
		std::vector<VTX_mindis> mindis;
		std::vector<VTX_base> base;
	};

	void read_vtxfile(std::string filename, std::vector<VTX>& vtx);
	//void write_vtxfile(std::string filename, std::vector<VTX>vtx);
	void write_vtxfile(std::string filepath, std::vector<VTX> vtx);

}


namespace netscan {

	class micro_track_t {
	public:
		double ax, ay;
		double x, y, z;
		double z1, z2;
		float px, py;
		int ph, pos, col, row, zone, isg;
		int64_t rawid;
	};

	class micro_track_subset_t {
	public:
		double ax, ay;
		double z;
		int ph;
		int pos, col, row, zone, isg;
		int64_t rawid;
	};

	class base_track_t {
	public:
		double ax, ay;
		double x, y, z;
		int pl;
		int isg, zone;
		int dmy;    // In ROOT, you will have to add this member because CINT does not handle 8byte alignment. 
		int64_t rawid;
		micro_track_subset_t m[2];
	};

	class linklet_t {
	public:
		int pos[2];
		double zproj;
		double dx, dy;
		double xc, yc;
		base_track_t b[2];

		// Below member is commented-out ( i.e. de-activated ) by default action.  
		// To access them, you need un-comment and give it to dump_linklet like --format 0 your-netscan-data-types-ui.h 
		//micro_track_t m[4];

	};

	bool read_linklet_txt(std::string filename, std::vector<linklet_t>& link, int output = 1);
	bool read_linklet_bin(std::string filename, std::vector<linklet_t>& link);
	bool read_basetrack_txt(std::string filename, std::vector<base_track_t>& base, int output = 1);
	void read_basetrack_extension(std::string filename, std::vector<base_track_t>& base, int PL, int zone);
	bool read_microtrack_txt(std::string filename, std::vector<micro_track_t>& micro, int output = 1);
	void read_microtrack_extension(std::string filename, std::vector<micro_track_t>& micro, int pos, int zone);

	void write_basetrack_txt(std::string filename, std::vector<base_track_t>& base);
	void write_basetrack_vxx(std::string filename, std::vector<base_track_t>& base, int PL, int zone);
	void write_basetrack_bin(std::string filename, std::vector<base_track_t>& base);
	void write_linklet_txt(std::string filename, std::vector<linklet_t> l);
	void write_linklet_bin(std::string filename, std::vector<linklet_t> l);
	void write_microtrack_txt(std::string filename, std::vector<micro_track_t>& micro);
	void write_microtrack_vxx(std::string filename, std::vector<micro_track_t>& micro, int pos, int zone);
	void write_microtrack_bin(std::string filename, std::vector<micro_track_t>& micro);


	base_track_t base_trans(base_track_t base, corrmap0::Corrmap corr);
	base_track_t base_trans_inv(base_track_t base, corrmap0::Corrmap corr);

}

//Chamber structure
namespace chamber0 {
	struct chamber_one_layer {
		double z;
		bool iron, water;
	};
	class chamber_structure {
	private:
		std::map<int, chamber_one_layer> chamber;
		/*
		pl1--iron--pl2--water--pl3--iron--pl4--water--pl5
			pl1 iron:true   water:false
			pl2 iron:false  water:ture
			pl3 iron:true   water:false
			pl4 iron:false  water:ture
		*/
	public:
		chamber_structure();
		chamber_structure(std::string);
		void disp_z_val();
		double dz(int pl1, int pl2);
		int number_of_iron(int pl1, int pl2);
		int number_of_water(int pl1, int pl2);
		double z_nominal(int pl);
	};

}
namespace chamber1 {
	class One_layer {
	public:
		double thick, radiation_length;
		std::string material;
	};
	class Chamber {
	public:
		std::map<int, double> z;
		std::vector<One_layer> layer;
	};
	void read_structure(std::string filneame, Chamber& chamber);
	std::map<int, double> base_z_convert(chamber1::Chamber& chamber);

}

namespace corrmap_3d {
	class align_param {
	public:
		int id, signal, ix, iy;
		//視野中心
		double x, y, z;
		//parameter(9);
		double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;

	};
	class align_param2 {
	public:
		align_param* corr_p[3];
		//3点の視野中心の重心(回転中心)
		double x, y, z;
		//parameter(9);
		double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;
	public:
		//3つのparameterから計算
		void Calc_9param();
		void Calc_9param_inv();


	};
	class align_param_3D {
	public:
		int id[3];
		double x[2][3], y[2][3], z[2][3];
		//parameter(9);
		double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;
	public:
		//3つのparameterから計算
		void Calc_9param();
		void Param_inv();

	};


	class Position {
	public:
		int id;
		double x, y, z;
	};

	bool triangle_internal_angle_cut(align_param2& param, double threshold_degree);

	std::vector<align_param> read_ali_param(std::string filename, bool output);
	std::map<int, std::vector<align_param>> read_ali_param_abs(std::string filename, bool output);

	std::vector <align_param2 >DelaunayDivide(std::vector <align_param >& corr);
	std::map<int, std::vector<align_param2>> DelaunayDivide_map(std::map<int, std::vector<align_param>>& corr);

	std::vector <std::pair<vxx::base_track_t*, align_param2*>>track_affineparam_correspondence(std::vector<vxx::base_track_t>& base, std::vector <align_param2>& param);
	std::vector <std::pair<vxx::base_track_t*, align_param2*>>track_affineparam_correspondence(std::vector<vxx::base_track_t*>& base, std::vector <align_param2>& param);
	align_param2* search_param(std::vector<align_param*>& param, vxx::base_track_t& base, std::multimap<int, align_param2*>& triangles);
	double select_triangle_vale(align_param2* param, vxx::base_track_t& base);
	void trans_base_all(std::vector < std::pair<vxx::base_track_t*, align_param2*>>& track_pair);
	void trans_base(std::vector<vxx::base_track_t*>& base, align_param2* param);
	void trans_basetrack(std::vector < vxx::base_track_t>& base, std::vector<align_param>& corr);
	void base_extra_PL0_to_PL1(std::vector<vxx::base_track_t>& base, std::vector<align_param2 >& param, bool coordinate_inverse);

	//param_3D
	std::vector <align_param_3D >Make_3D_Param(std::vector <align_param >& corr, double inner_angle_thr = 20);
	//basetrack-alignment mapの対応
	std::vector <std::pair<vxx::base_track_t*, align_param_3D*>>track_affineparam_correspondence(std::vector<vxx::base_track_t>& base, std::vector <align_param_3D>& param);
	align_param_3D* search_param(std::vector<Position>& param, vxx::base_track_t& base, std::multimap<int, align_param_3D*>& triangles);
	double select_triangle_vale(align_param_3D* param, vxx::base_track_t& base);
	//変換 zshrink補正-->9para変換
	void trans_base_all(std::vector < std::pair<vxx::base_track_t*, align_param_3D*>>& track_pair);
	void trans_base(std::vector<vxx::base_track_t*>& base, align_param_3D* param);

	//basetrack-alignment mapの対応
	std::vector <std::pair<mfile0::M_Base*, align_param2*>>track_affineparam_correspondence(std::vector<mfile0::M_Base*>& base, std::vector <align_param2>& param);
	align_param2* search_param(std::vector<Position>& param, mfile0::M_Base* base, std::multimap<int, align_param2*>& triangles);
	double select_triangle_vale(align_param2* param, mfile0::M_Base* base);
	//変換 zshrink補正-->9para変換
	void trans_base_all(std::vector < std::pair<mfile0::M_Base*, align_param2*>>& track_pair);
	void trans_base(std::vector<mfile0::M_Base*>& base, align_param2* param);


}

namespace chain {
	class Header {
	public:
		int n_pos;
		std::vector<int> use_pl;
		int PmPeke;
		//#       N_pos = 130
		//	#       UsePos = 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690 700 710 720 730 740 750 760 770 780 790 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 1110 1120 1130 1140 1150 1160 1170 1180 1190 1200 1210 1220 1230 1240 1250 1260 1270 1280 1290 1300 1310 1320 1330
		//	#       PmPeke = 5
	};
	class Group_header {
	public:
		int64_t group_id, nchain;
	};
	class Group_header_all :public Group_header {
	public:
		int64_t nseg;
		int start_pl, end_pl, root_num, leaf_num;
		double nseg_ave, nseg_sig;
		//$1 = Group ID
		//	$2 = Groupに属するchain数
		//	$3 = Groupのstart plate（最先端segの、UsePosで設定した最初のplate番号を1としたplate番号。Pos / 10 ではない）
		//	$4 = Groupのend plate（最終端segの、UsePosで設定した最初のplate番号を1としたplate番号。Pos / 10 ではない）
		//	$5 = Groupが持つsegment数
		//	$6 = rootの数（先端のseg数）
		//	$7 = leafの数（終端のseg数）
		//	$8 = 1plate内のsegment数の平均値（ペケもsegmentとしてカウント）
		//	$9 = $8のσ

	};
	class Chain {
	public:
		int64_t chain_id;
		int start_pl, end_pl, nseg;
		std::vector<std::pair<int, int>> rawid;
	};
	class Chain_all :public Chain {
	public:
		int root_leaf_id, root_id, start_pl;
		std::vector<int> peke_num;
		//chainのheader（段が1段下がっている行）
		//	$1 = chain ID
		//	$2 = root と leaf が同じ chain に対して同じ値になる ID
		//	$3 = root が同じ chain に対して同じ値になる ID
		//	$4 = start plate(rootのplate)
		//	$5 = end plate（leafのplate）
		//	$6 = nseg
		//	$7 = 1ペケの数　　　　許すペケの数 PmPeke = 1 以上の指定で出力される
		//	$8 = 2ペケの数        許すペケの数 PmPeke = 2 以上の指定で出力される
	};
	class Footer {
	public:
		int64_t n_linklet, n_folded_linklet, n_chain, n_gcomm, n_groot, n_gends;
		//#    N_linklet = 2597
		//	#    N_folded - linklet = 434
		//	#    N_chain = 1788
		//	#    N_gcomm = 1576
		//	#    N_groot = 1664
		//	#    N_gends = 1774
	};
	class Chain_file {
	public:
		Header header;
		Footer footer;
		std::vector<Group_header> groups;
		std::vector<std::vector<Chain>> chains;
	};

	void read_chain(std::string file_path, Chain_file& chain);
}

namespace Sharing_file {
	class Sharing_file {
	public:
		int32_t pl, ecc_id, oss_id, fixedwall_id, trackerwall_id, spotid, zone[2], rawid[2], unix_time, tracker_track_id, babymind_bunch, babymind_nplane, charge, entry_in_daily_file, eventid, track_type, ecc_track_type;
		//spotid:spotA * 100 + spotB
		float chi2_shifter[4], babymind_momentum;
		//spotid:spotA * 100 + spotB
		// chi2_shifter : [0]:ECC-fixedwall [1]:fixedwall-TSS [2]:TSS-tracker
		static bool sort_unix_time(const Sharing_file& lhs, const Sharing_file& rhs) {
			if (lhs.unix_time == rhs.unix_time)return lhs.tracker_track_id < rhs.tracker_track_id;
			return lhs.unix_time < rhs.unix_time;
		}
		static bool sort_eventid(const Sharing_file& lhs, const Sharing_file& rhs) {
			if (lhs.eventid == rhs.eventid)return sort_unix_time(lhs, rhs);
			//lhs.unix_time < rhs.unix_time;
			return lhs.eventid < rhs.eventid;
		}
		std::string Print_content();
	};

	std::vector<Sharing_file> Read_sharing_file_bin(std::string filename);
	std::vector<Sharing_file> Read_sharing_file_txt(std::string filename);
	void Write_sharing_file_bin(std::string filename, std::vector<Sharing_file>& sharing_file_v);
	void Write_sharing_file_txt(std::string filename, std::vector<Sharing_file>& sharing_file_v);
	std::vector<Sharing_file> Read_sharing_file_extension(std::string filename);

}
namespace Momentum_recon {

	class microtrack_minimum {
	public:
		int ph, zone, view, imager, pixelnum, hitnum;
		microtrack_minimum();
	};
	class Mom_basetrack {
	public:
		int pl, rawid;
		float ax, ay, x, y, z;
		microtrack_minimum m[2];
		Mom_basetrack();
	};
	class Mom_chain {
	public:
		//stop_flg=0 : penetarte / 1:bm stop / 2:ecc stop
		int chainid, stop_flg, particle_flg, direction, charge_sign;
		//[0]:muon [1]:proton
		double ecc_range_mom[2], ecc_mcs_mom[2], bm_range_mom, bm_curvature_mom;
		double ecc_range_mom_error[2][2], ecc_mcs_mom_error[2][2], bm_range_mom_error[2], bm_curvature_mom_error[2];
		//muon_likelihood   =(vph-  muon_mean)/  muon_sigma (符号付き)
		//proton_likelihood =(vph-proton_mean)/proton_sigma (符号付き)
		double muon_likelihood, proton_likelihood;
		//mfile座標系
		std::vector<Mom_basetrack> base;
		//localな座標系
		std::vector<std::pair<Mom_basetrack, Mom_basetrack>> base_pair;

		//コンストラクタ
		Mom_chain();
		double Get_muon_mcs_pb();
		void Get_muon_pb_mcs_error(double* pb_error);
		double Get_proton_mcs_pb();
		void Get_proton_pb_mcs_error(double* pb_error);
	};
	class Event_information {
	public:
		int groupid, unix_time, tracker_track_id, entry_in_daily_file, vertex_pl, ECC_id, vertex_material;
		double vertex_position[3], true_vertex_position[3], weight, nu_energy, nu_ax, nu_ay;
		std::vector<Mom_chain> chains;
		std::vector<Mom_chain> true_chains;

		//コンストラクタ
		Event_information();
		//0:water,2:iron -1:ini
	};

	void Write_Event_information_header(std::ofstream& ofs, Event_information& ev);
	bool Read_Event_information_header(std::ifstream& ifs, Event_information& ev, int& chian_num, int& true_chian_num);
	void Write_mom_chain_header(std::ofstream& ofs, Mom_chain& mom_chain);
	bool Read_mom_chain_header(std::ifstream& ifs, Mom_chain& mom_chain, int& base_num, int& base_pair_num);
	void Write_Event_information_bin(std::string filename, std::vector<Event_information>& ev_v);
	std::vector<Event_information> Read_Event_information_bin(std::string filename);
	void Write_Event_information_bin_block(std::ofstream& ofs, std::vector<Event_information>& ev_v);
	std::vector<Event_information> Read_Event_information_bin_block(std::ifstream& ifs, int max_event_num);

	void Write_Event_information_txt(std::string filename, std::vector<Event_information>& ev_v);
	std::vector<Event_information> Read_Event_information_txt(std::string filename);
	void input_basetrack_information(std::vector<Mom_basetrack>& base, std::vector<std::pair<Mom_basetrack, Mom_basetrack>>& base_pair);
	std::vector<Event_information> Read_Event_information_extension(std::string filename);
	void Write_Event_information_extension(std::string filename, std::vector<Event_information>& ev_v);
}
