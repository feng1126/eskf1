/*
 * @Author: your name
 * @Date: 2021-07-26 15:16:42
 * @LastEditTime: 2021-07-28 16:13:49
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \code\include\fusionalgorithm\enu.h
 */
#pragma once
#include <math.h>

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#define PI 3.141592653589793
class ENU
{
public:
	ENU();
	void setStart(const double& a, const double& b, const float& c)
	{
		m_lla_ref[0] = a / 180.0 * PI;
		m_lla_ref[1] = b / 180.0 * PI;
		m_lla_ref[2] = c;
		pos2ecef(m_lla_ref, m_ecef_ref);
	}

	void  getENU(const double& a, const double& b, const float& c, double& i,  double& j, double& k)
	{
		double m[3],n[3], r[3],ecef[3] = { 0 } ;
		m[0] = a / 180.0 * PI;
		m[1] = b / 180.0 * PI;
		m[2] = c;
		pos2ecef(m, ecef);
		r[0] = ecef[0] - m_ecef_ref[0];
		r[1] = ecef[1] - m_ecef_ref[1];
		r[2] = ecef[2] - m_ecef_ref[2];
		ecef2enu(m_lla_ref, r, n);
		i = n[0];
		j = n[1];
		k = n[2];
		return ;
	}

	void getLLA(const double& a, const double& b, const float& c ,double& i, double& j,double& k )
	{
		double n[3], r[3], lla[3];
		r[0] = a;
		r[1] = b;
		r[2] = c;

		enu2ecef(m_lla_ref, r, n);
		n[0] = n[0] + m_ecef_ref[0];
		n[1] = n[1] + m_ecef_ref[1];
		n[2] = n[2] + m_ecef_ref[2];
		ecef2pos(n, lla);
		i = lla[0] * 180 / PI;
		j = lla[1] * 180 / PI;
		k = lla[2];
		return ;
	}

private:
	double m_ecef_ref[3], m_lla_ref[3];
	double dot(const double* a, const double* b, int n);
	void ecef2pos(const double* r, double* pos);
	void pos2ecef(const double* pos, double* r);
	void xyz2enu(const double* pos, double* E);
	void ecef2enu(const double* pos, const double* r, double* e);
	void enu2ecef(const double* pos, const double* e, double* r);
	void matmul(const char* tr, int n, int k, int m, double alpha,const double* A, const double* B, double beta, double* C);
};
