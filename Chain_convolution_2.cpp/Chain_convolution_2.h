#pragma once
#include <FILE_structure.h>


class Chain_path {
	std::vector<mfile0::M_Base> base;
	double ax, ay, x, y, z;
	double length_straight, length_path;

public:
	Chain_path();
	Chain_path(mfile0::M_Base &b0, mfile0::M_Base&b1);
	Chain_path(std::vector<mfile0::M_Base>&b);
	Chain_path(std::set<mfile0::M_Base>&b);

	bool operator==(const Chain_path& rhs)const;
	bool operator<(const Chain_path& rhs)const;
	bool operator>(const Chain_path& rhs)const;


	void Line_Fit();
	void Line_Fit(int pl0, int pl1);
	bool Add_Basetrack(mfile0::M_Base &b);
	void Calc_length();
	std::pair<mfile0::M_Base, mfile0::M_Base> Get_path_edge();
	std::vector<mfile0::M_Base> Get_all_base();

	double Get_line_ax() { return ax; };
	double Get_line_ay() { return ay; };
	double Get_line_x(double z_pos) { return ax * (z_pos - z) + x; };
	double Get_line_y(double z_pos) { return ay * (z_pos - z) + y; };
};
bool sort_chain_base(const mfile0::M_Base &left, mfile0::M_Base &right) {
	if (left.pos == right.pos)return left.rawid < right.rawid;
	return left.pos < right.pos;
}
Chain_path::Chain_path() {
	base.clear();
	ax = 0; 
	ay = 0; 
	x = 0; 
	y = 0; 
	z = 0;
	length_path = 0;
	length_straight = 0;
}

Chain_path::Chain_path(mfile0::M_Base &b0, mfile0::M_Base&b1) {
	Add_Basetrack(b0);
	Add_Basetrack(b1);
	Calc_length();
}
Chain_path::Chain_path(std::vector<mfile0::M_Base>&b) {
	for (int i = 0; i < b.size(); i++) {
		Add_Basetrack(b[i]);
	}
	Calc_length();
}
Chain_path::Chain_path(std::set<mfile0::M_Base>&b) {
	for (auto itr = b.begin(); itr != b.end();itr++) {
		mfile0::M_Base b_in = *itr;
		Add_Basetrack(b_in);
	}
	Calc_length();
}

bool Chain_path::operator==(const Chain_path& rhs)const {
	int pos0[2] = { this->base.begin()->pos, this->base.rbegin()->pos };
	int pos1[2] = { rhs.base.begin()->pos, rhs.base.rbegin()->pos };
	int raw0[2] = { this->base.begin()->rawid, this->base.rbegin()->rawid };
	int raw1[2] = { rhs.base.begin()->rawid, rhs.base.rbegin()->rawid };
	//printf("operator ==\n");
	for (int i = 0; i < 2; i++) {
		if (pos0[i] != pos1[i])return false;
		if (raw0[i] != raw1[i])return false;
	}
	return true;
}
bool Chain_path::operator<(const Chain_path& rhs) const {
	int pos0[2] = { this->base.begin()->pos, this->base.rbegin()->pos };
	int pos1[2] = { rhs.base.begin()->pos, rhs.base.rbegin()->pos };
	int raw0[2] = { this->base.begin()->rawid, this->base.rbegin()->rawid };
	int raw1[2] = { rhs.base.begin()->rawid, rhs.base.rbegin()->rawid };
	bool flg = false;
	for (int i = 0; i < 2; i++) {
		if (pos0[i] < pos1[i])return true;
		else if (pos0[i] > pos1[i])return false;
	}
	for (int i = 0; i < 2; i++) {
		if (raw0[i] < raw1[i])return true;
		else if (raw0[i] > raw1[i])return false;
	}
	return false;

	if (flg) {
		printf("(%d,%d,%d,%d)<(%d,%d,%d,%d)\n", pos0[0], pos0[1], raw0[0], raw0[1], pos1[0], pos1[1], raw1[0], raw1[1]);
	}
	else {
		printf("(%d,%d,%d,%d)>(%d,%d,%d,%d)\n", pos0[0], pos0[1], raw0[0], raw0[1], pos1[0], pos1[1], raw1[0], raw1[1]);
	}
}

bool Chain_path::operator>(const Chain_path& rhs) const {
	int pos0[2] = { base.begin()->pos, base.rbegin()->pos };
	int pos1[2] = { rhs.base.begin()->pos, rhs.base.rbegin()->pos };
	int raw0[2] = { base.begin()->rawid, base.rbegin()->rawid };
	int raw1[2] = { rhs.base.begin()->rawid, rhs.base.rbegin()->rawid };
	for (int i = 0; i < 2; i++) {
		if (pos0[i] > pos1[i])return true;
	}
	for (int i = 0; i < 2; i++) {
		if (raw0[i] > raw1[i])return true;
	}
	return false;
}


bool Chain_path::Add_Basetrack(mfile0::M_Base &b) {
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		if (itr->pos == b.pos)return false;
	}
	base.push_back(b);
	sort(base.begin(), base.end(),sort_chain_base);
}
void Chain_path::Line_Fit() {
	Line_Fit(base.begin()->pos / 10, base.rbegin()->pos / 10);
}

void Chain_path::Line_Fit(int pl0,int pl1) {
	double sum_x, sum_y, sum_z, sum_zz, sum_xz, sum_yz, n;
	sum_x = 0;
	sum_y = 0;
	sum_z = 0;
	sum_zz = 0;
	sum_xz = 0;
	sum_yz = 0;
	n = 0;
	int pl;
	for (int i = 0; i < base.size(); i++) {
		pl = base[i].pos / 10;
		if (pl < pl0)continue;
		if (pl1 < pl)continue;
		sum_x += base[i].x;
		sum_y += base[i].y;
		sum_z += base[i].z;
		sum_zz += base[i].z*base[i].z;
		sum_xz += base[i].x*base[i].z;
		sum_yz += base[i].y*base[i].z;
		n += 1;
	}

	double denominator;
	denominator = n * sum_zz - sum_z * sum_z;
	if (n ==0) {
		x = 0;
		y = 0;
		z = 0;
		ax = 0;
		ay = 0;
	}
	else if(n==1){
		for (int i = 0; i < base.size(); i++) {
			pl = base[i].pos / 10;
			if (pl < pl0)continue;
			if (pl1 < pl)continue;
			x = base[i].x;
			y = base[i].y;
			z = base[i].z;
			ax = base[i].ax;
			ay = base[i].ay;
		}
	}
	else if (n == 2) {
		x = sum_x / 2;
		y = sum_y / 2;
		z = sum_z / 2;
		double dx, dy, dz;
		n = 0;
		for (int i = 0; i < base.size(); i++) {
			pl = base[i].pos / 10;
			if (pl < pl0)continue;
			if (pl1 < pl)continue;
			if (n == 0) {
				dx = base[i].x;
				dy = base[i].y;
				dz = base[i].z;
			}
			else if (n == 1) {
				dx -= base[i].x;
				dy -= base[i].y;
				dz -= base[i].z;
			}
			n++;
		}
		if (dz < 0.000001) {
			ax = 9999999;
			ay=999999;
		}
		else {
			ax = dx / dz;
			ay = dy / dz;
		}
	}
	else if (denominator < 0.001) {
		double dx, dy, dz;
		for (int i = 0; i < base.size(); i++) {
			pl = base[i].pos / 10;
			if (pl >= pl0) {
				x = base[i].x;
				y = base[i].y;
				z = base[i].z;

				dx = base[i].x;
				dy = base[i].y;
				dz = base[i].z;
				break;
			}
		}
		for (int i = base.size()-1; i>=0; i--) {
			pl = base[i].pos / 10;
			if (pl <= pl1) {
				x += base[i].x;
				y += base[i].y;
				z += base[i].z;

				dx -= base[i].x;
				dy -= base[i].y;
				dz -= base[i].z;
				break;
			}

		}
		x /= 2;
		y /= 2;
		z /= 2;
		if (dz < 0.000001) {
			ax = 9999999;
			ay = 999999;
		}
		else {
			ax = dx / dz;
			ay = dy / dz;
		}

	}
	else {
		ax = (n*sum_xz - sum_x * sum_z) / denominator;
		ay = (n*sum_yz - sum_y * sum_z) / denominator;
		x = (sum_zz*sum_x - sum_xz * sum_z) / denominator;
		y = (sum_zz*sum_y - sum_yz * sum_z) / denominator;
		z = 0;
	}

}
void Chain_path::Calc_length() {
	length_path = 0;

	for (int i = 0; i < base.size();i++) {
		if (i + 1 == base.size())continue;
		length_path += sqrt(pow(base[i].x - base[i + 1].x, 2) + pow(base[i].y - base[i + 1].y, 2) + pow(base[i].z - base[i + 1].z, 2));
	}
	Line_Fit();
	double z0 = base.begin()->z;
	double z1 = base.rbegin()->z;
	double dx, dy, dz;
	dx = ax * (z0 - z1);
	dy = ay * (z0 - z1);
	dz =  (z0 - z1);
	length_straight = sqrt(dx*dx + dy * dy + dz * dz);
}
std::pair<mfile0::M_Base, mfile0::M_Base> Chain_path::Get_path_edge() {
	return std::make_pair(base[0], base[base.size() - 1]);
}

std::vector<mfile0::M_Base> Chain_path::Get_all_base() {
	return base;
}

bool operator==(const std::pair<Chain_path, Chain_path>&left, const std::pair<Chain_path, Chain_path>&right) {
	if (left.first == right.first && left.second == right.second)return true;
	return false;
}

void Print_path(mfile0::M_Base &b0, mfile0::M_Base &b1) {
	printf("(%d,%d)-(%d,%d)", b0.pos, b0.rawid, b1.pos, b1.rawid);
}void Print_path(const mfile0::M_Base &b0,const mfile0::M_Base &b1) {
	printf("(%d,%d)-(%d,%d)", b0.pos, b0.rawid, b1.pos, b1.rawid);
}