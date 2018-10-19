#include "BezierCurve.h"
#include "NaturalSpline.h"
#include "Logger.h"
#include <gtest/gtest.h>
#include <random>
#include <fstream>

TEST(bezier_spline_test,BezierSpline)
{
	BezierPolyLine<float> p;

	p.Append({0,0},false);
	p.Append({200,50},false);
	p.Insert(1,{150,100},false);
	p.Insert(1,{100,60},false);
	p.Insert(1,{50,100},false);
	for (std::size_t k=0; k<p.GetNumNodes()-1; ++k)
		p.MakeNodeSmooth(k);//p.MakeStraight(k);
	//p.Dump();

	/*
	ConstraintSet<2,double> CS {
		{ {{1,2}} , 3 },
		{ {{3,4}} , 1 }
	};
	printf("%g\n", CS(0,1));
	return 0;
	*/

	//p.MoveLeftCtrlTo(2,60,50);
	//p.MoveRightCtrlTo(2,140,50);
	//p.MakeSmooth(0);
	//p.MakeSmooth(1);
	/*p.MakeSmooth(2);
	p.MakeSmooth(3);
	p.MakeSmooth(4);
	p.MakeStraight(0);*/

	BezierCurve<float> S(p);
	printf("number of segments: %u\n",S.GetNumSegments());
	float t_min = S.ProjectOn(std::array<float,2>{{100.0f,40.0f}});
	auto p_prj = S.Eval(t_min);
	printf("t_min = %g => (%g,%g)\n", t_min, p_prj[0], p_prj[1]);
//return;

	for (float t : Lerp<float>(S.GetMinParam(),S.GetMaxParam(),100))
	{
		std::array<float,2> v = S.Eval(t);
		for (float s : v)
			printf(" %g", s);
		printf("\n");
	}
}

TEST(bezier_spline,Copy)
{
	// TODO
}

TEST(bezier_spline_test,Conversion)
{
	// Spline
	double X[] = { 0, 1, 2, 3, 4 };
	double Y[] = { -10, 5, 3, -1, 5};
	NaturalSpline<double> spline =
		NaturalSpline<double>::Interpolate(X,Y,5,3);

	// Konvertierung des Spline in eine Poly-Line
	BezierPolyLine<double> bpoly(spline);
	bpoly.Dump();

	// Konvertierung der Poly-Line in eine Bèzier-Kurve
	BezierCurve<double> bcurve(bpoly);

	printf("parameter interval t = [%g,%g]\n",
		bcurve.GetMinParam(), bcurve.GetMaxParam());

	// y-Komponente der Bèzier-Kurve in einen
	//const Curve<double>& spline_y = bcurve[1];

	std::ofstream ofs("bezier.csv");
	for (double t : Lerp<double>(bcurve.GetMinParam(),bcurve.GetMaxParam(),100))
	{
		std::array<double,2> p = bcurve.Eval(t);
		double y = spline.Eval(t);
		ofs << p[0] << ' ' << p[1] << ' ' << y << std::endl;
		EXPECT_LT(std::abs(p[1]-y),1e-6);
	}

	printf("reconstructed poly line (y(x))\n");
	BezierPolyLine<double> bpoly2(bcurve[1]);
	bpoly2.Dump();

	printf("reconstructed poly line (x(t),y(t))\n");
	BezierPolyLine<double> bpoly3(bcurve);
	bpoly3.Dump();
}

TEST(bezier_spline_test,Conversion2)
{
	/*
	BezierPolyLine<double> bpoly;

	bpoly.Append(100,30,true);
	bpoly.Append(150,-10,true);

	bpoly.MoveRightCtrlTo(0,120,30);
	bpoly.MoveLeftCtrlTo(1,130,-10);

	BezierCurve<double> bspl(bpoly);

	std::ofstream ofs("conv2.csv");
	for (double t : Lerp<double>(bspl.GetMinParam(),bspl.GetMaxParam(),100))
	{
		auto p = bspl.Eval(t);
		ofs << p[0] << ' ' << p[1] << std::endl;
	}

	Curve<double>& crv_x = bspl[0];
	Curve<double>& crv_y = bspl[1];

	for (std::size_t s=0; s<crv_x.GetNumSegments(); ++s)
	{
		printf("s=%zi: X: c0=%g, c1=%g, c2=%g, c3=%g\n",
			s, crv_x.GetCoeff(s,0), crv_x.GetCoeff(s,1), crv_x.GetCoeff(s,2), crv_x.GetCoeff(s,3));
		printf("s=%zi: Y: c0=%g, c1=%g, c2=%g, c3=%g\n",
			s, crv_y.GetCoeff(s,0), crv_y.GetCoeff(s,1), crv_y.GetCoeff(s,2), crv_y.GetCoeff(s,3));
	}
	*/
#if 1
	// Spline
	double X[] = { -10, 5, 13, 47, 50 };
	//double X[] = { 0, 1, 2, 3 };
	double Y[] = { -10, -50, 13, -60, 0 };
	NaturalSpline<double> spline =
		NaturalSpline<double>::Interpolate(X,Y,5,3);

	// Konvertierung des Spline in eine Poly-Line
	BezierPolyLine<double> bpoly(spline);
	bpoly.Dump();

	std::ofstream ofs("bezier2.csv");
	//std::ofstream ofss("bezier3.csv");
	BezierCurve<double> bezier(bpoly);
	for (double t : Lerp<double>(bezier.GetMinParam(),bezier.GetMaxParam(),100))
	{
		auto p = bezier.Eval(t);
		ofs << p[0] << ' ' << p[1] << std::endl;
		//ofss << t << ' ' << spline(60*t) << std::endl;
		//printf("%g %g\n",p[0],p[1]);
	}
#endif
}

TEST(bezier_spline_test,Subdivision)
{
	// Spline
	double X[] = { 0, 1, 2, 3, 4 };
	double Y[] = { -10, 5, 3, -1, 5};
	NaturalSpline<double> spline =
		NaturalSpline<double>::Interpolate(X,Y,5,3);

	// Konvertierung des Spline in eine Poly-Line
	BezierPolyLine<double> bpoly(spline);
	bpoly.Dump();

	// Konvertierung der Poly-Line in eine Bèzier-Kurve
	BezierCurve<double> bcurve(bpoly);

	printf("parameter interval t = [%g,%g]\n",
		bcurve.GetMinParam(), bcurve.GetMaxParam());

	bpoly.Subdivide(0,0.7);
	BezierCurve<double> bcurve_subdiv(bpoly);

	printf("parameter interval t = [%g,%g]\n",
		bcurve_subdiv.GetMinParam(), bcurve_subdiv.GetMaxParam());

	bpoly.Dump();

	std::ofstream ofs("subdiv.csv");
	LinearMap<double> mp_bcurve(0,99,bcurve.GetMinParam(),bcurve.GetMaxParam());
	LinearMap<double> mp_bcurve_subdiv(0,99,bcurve_subdiv.GetMinParam(),bcurve_subdiv.GetMaxParam());
	for (int ti=0; ti<100; ++ti)
	{
		std::array<double,2> p = bcurve.Eval(mp_bcurve(ti));
		std::array<double,2> q = bcurve_subdiv.Eval(mp_bcurve_subdiv(ti));

		double t_prj = bcurve_subdiv.ProjectOn(p);
		//printf("t_prj=%g vs. %g, d=%g\n",
		//	t_prj, mp_bcurve_subdiv(ti), std::abs(t_prj-mp_bcurve_subdiv(ti)));
		//EXPECT_LT(std::abs(p[1]-y),1e-6);

		std::array<double,2> r = bcurve_subdiv.Eval(t_prj);

		ofs << p[0] << ' ' << p[1] << ' ' << q[0] << ' ' << q[1]
			<< ' ' << r[0] << ' ' << r[1] << std::endl;
	}
	ofs.close();

	std::ofstream nfs("nodes.csv");
	for (std::size_t k=0; k<bpoly.GetNumNodes(); ++k)
	{
		nfs << bpoly.Get(k)[0] << ' ' << bpoly.Get(k)[1] << std::endl;
	}
}

TEST(bezier_spline_test,MoveSegment)
{
	BezierPolyLine<double> bpoly;

	bpoly.Append({85,241},false);
	bpoly.Append({210,160},false);
	bpoly.MoveRightCtrlTo(0,{20,155});
	bpoly.MoveLeftCtrlTo(1,{193,85});

	bpoly.Dump();
	bpoly.MoveSegmentTo(0.8,{150,120});
	bpoly.Dump();

	BezierCurve<double> bcurve(bpoly);
	Vector<double,2> v = bcurve.Eval(0.8);

	printf("(%g,%g)\n", v[0],v[1]);
}

TEST(bezier_spline_test,Render)
{
	// Spline
	double X[] = { 0, 1, 2, 3, 4 };
	double Y[] = { -10, 5, 3, -1, 5};
	NaturalSpline<double> spline =
		NaturalSpline<double>::Interpolate(X,Y,5,3);

	// Konvertierung des Spline in eine Poly-Line
	BezierPolyLine<double> bpoly(spline);

	std::vector<Vector<double,2>> path = bpoly.ToPath(0.01);

	std::ofstream ofs("rendering.csv");
	for (Vector<double,2>& v : path)
	{
		ofs << v[0] << ' ' << v[1] << std::endl;
	}
	ofs.close();
}

