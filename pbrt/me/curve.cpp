
// shapes/curve.cpp*
#include "me/curve.h"
#include "paramset.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Curves", curveBytes);
STAT_PERCENT("Intersections/Ray-curve intersection tests", nHits, nTests);
STAT_INT_DISTRIBUTION("Intersections/Curve refinement level", refinementLevel);
STAT_COUNTER("Scene/Curves", nCurves);
STAT_COUNTER("Scene/Split curves", nSplitCurves);

// Curve Utility Functions
static Point3f BlossomBezier(const Point3f p[4], Float u0, Float u1, Float u2) {
	Point3f a[3] = {
		Lerp(u0, p[0], p[1]),
		Lerp(u0, p[1], p[2]),
		Lerp(u0, p[2], p[3])
	};
	Point3f b[2] = {
		Lerp(u1, a[0], a[1]),
		Lerp(u1, a[1], a[2])
	};
	
	return Lerp(u2, b[0], b[1]);
}

inline void SubdivideBezier(const Point3f cp[4], Point3f cpSplit[7]) {
	cpSplit[0] = cp[0];
	cpSplit[1] = (cp[0] + cp[1]) / 2;
	cpSplit[2] = (cp[0] + 2 * cp[1] + cp[2]) / 4;
	cpSplit[3] = (cp[0] + 3 * cp[1] + 3 * cp[2] + cp[3]) / 8;
	cpSplit[4] = (cp[1] + 2 * cp[2] + cp[3]) / 4;
	cpSplit[5] = (cp[2] + cp[3]) / 2;
	cpSplit[6] = cp[3];
}

static Point3f EvalBezier(const Point3f cp[4], Float u, Vector3f* deriv = nullptr) {
	Point3f cp1[3] = {
		Lerp(u, cp[0], cp[1]),
		Lerp(u, cp[1], cp[2]),
		Lerp(u, cp[2], cp[3])
	};
	Point3f cp2[2] = {
		Lerp(u, cp1[0], cp1[1]),
		Lerp(u, cp1[1], cp1[2])
	};

	if (deriv) {
		if ((cp2[1] - cp2[0]).LengthSquared() > 0) {
			*deriv = 3 * (cp2[1] - cp2[0]);
		} else {
			*deriv  = cp[3] - cp[0];
		}
	}
	
	return Lerp(u, cp2[0], cp2[1]);
}



CurveCommon::CurveCommon(const Point3f c[4], Float w0, Float w1, CurveType type, const Normal3f* norm)
	: type(type){
	width[0] = w0;
	width[1] = w1;
	for (int i = 0; i < 4; ++i) {
		cpObj[i] = c[i];
	}
	if (norm) {
		n[0] = Normalize(norm[0]);
		n[1] = Normalize(norm[1]);
		normalAngle = std::acos(Clamp(Dot(n[0], n[1]), 0, 1));
		invSinNormalAngle = 1 / std::sin(normalAngle);
	}
	++nCurves;
}

//for one curvecom, segments
std::vector<std::shared_ptr<Shape>> CreateCurve(
	const Transform* o2w, const Transform* w2o, bool reverseOrientation,
	const Point3f* c, Float w0, Float w1, CurveType type,
	const Normal3f* norm, int splitDepth) {
	std::vector<std::shared_ptr<Shape>> segments;
	std::shared_ptr<CurveCommon> common = std::make_shared<CurveCommon>(c, w0, w1, type, norm);
	const int nSegments = 1 << splitDepth;
	const Float fN = static_cast<Float>(nSegments);
	segments.reserve(nSegments);
	for (int i = 0; i < nSegments; ++i) {
		Float uMin = i / fN;
		Float uMax = (i + 1) / fN;
		segments.push_back(std::make_shared<Curve>(o2w, w2o, reverseOrientation, common, uMin,uMax));
		++nSplitCurves;
	}
	curveBytes += sizeof(CurveCommon) + nSegments * sizeof(Curve);
	return segments;
}





Interaction Curve::Sample(const Point2f &u, Float *pdf) const {
    LOG(FATAL) << "Curve::Sample not implemented.";
    return Interaction();
}

std::vector<std::shared_ptr<Shape>> CreateCurveShape(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params) {
    Float width = params.FindOneFloat("width", 1.f);
    Float width0 = params.FindOneFloat("width0", width);
    Float width1 = params.FindOneFloat("width1", width);

    int degree = params.FindOneInt("degree", 3);
    if (degree != 2 && degree != 3) {
        Error("Invalid degree %d: only degree 2 and 3 curves are supported.",
              degree);
        return {};
    }

    std::string basis = params.FindOneString("basis", "bezier");
    if (basis != "bezier" && basis != "bspline") {
        Error("Invalid basis \"%s\": only \"bezier\" and \"bspline\" are "
              "supported.", basis.c_str());
        return {};
    }

    int ncp;
    const Point3f *cp = params.FindPoint3f("P", &ncp);
    int nSegments;
    if (basis == "bezier") {
        // After the first segment, which uses degree+1 control points,
        // subsequent segments reuse the last control point of the previous
        // one and then use degree more control points.
        if (((ncp - 1 - degree) % degree) != 0) {
            Error("Invalid number of control points %d: for the degree %d "
                  "Bezier basis %d + n * %d are required, for n >= 0.", ncp,
                  degree, degree + 1, degree);
            return {};
        }
        nSegments = (ncp - 1) / degree;
    } else {
        if (ncp < degree + 1) {
            Error("Invalid number of control points %d: for the degree %d "
                  "b-spline basis, must have >= %d.", ncp, degree, degree + 1);
            return {};
        }
        nSegments = ncp - degree;
    }


    CurveType type;
    std::string curveType = params.FindOneString("type", "flat");
    if (curveType == "flat")
        type = CurveType::Flat;
    else if (curveType == "ribbon")
        type = CurveType::Ribbon;
    else if (curveType == "cylinder")
        type = CurveType::Cylinder;
    else {
        Error("Unknown curve type \"%s\".  Using \"cylinder\".", curveType.c_str());
        type = CurveType::Cylinder;
    }

    int nnorm;
    const Normal3f *n = params.FindNormal3f("N", &nnorm);
    if (n != nullptr) {
        if (type != CurveType::Ribbon) {
            Warning("Curve normals are only used with \"ribbon\" type curves.");
            n = nullptr;
        } else if (nnorm != nSegments + 1) {
            Error(
                "Invalid number of normals %d: must provide %d normals for ribbon "
                "curves with %d segments.", nnorm, nSegments + 1, nSegments);
            return {};
        }
    } else if (type == CurveType::Ribbon) {
        Error(
            "Must provide normals \"N\" at curve endpoints with ribbon "
            "curves.");
        return {};
    }

    int sd = params.FindOneInt("splitdepth",
                               int(params.FindOneFloat("splitdepth", 3)));

    std::vector<std::shared_ptr<Shape>> curves;
    // Pointer to the first control point for the current segment. This is
    // updated after each loop iteration depending on the current basis.
    const Point3f *cpBase = cp;
	//for each curvecommon, segmeng 
    for (int seg = 0; seg < nSegments; ++seg) {
        Point3f segCpBezier[4];

        // First, compute the cubic Bezier control points for the current
        // segment and store them in segCpBezier. (It is admittedly
        // wasteful storage-wise to turn b-splines into Bezier segments and
        // wasteful computationally to turn quadratic curves into cubics,
        // but yolo.)
        if (basis == "bezier") {
            if (degree == 2) {
                // Elevate to degree 3.
                segCpBezier[0] = cpBase[0];
                segCpBezier[1] = Lerp(2.f/3.f, cpBase[0], cpBase[1]);
                segCpBezier[2] = Lerp(1.f/3.f, cpBase[1], cpBase[2]);
                segCpBezier[3] = cpBase[2];
            } else {
                // Allset.
                for (int i = 0; i < 4; ++i)
                    segCpBezier[i] = cpBase[i];
            }
            cpBase += degree;
        } else {
            // Uniform b-spline.
            if (degree == 2) {
                // First compute equivalent Bezier control points via some
                // blossiming.  We have three control points and a uniform
                // knot vector; we'll label the points p01, p12, and p23.
                // We want the Bezier control points of the equivalent
                // curve, which are p11, p12, and p22.
                Point3f p01 = cpBase[0];
                Point3f p12 = cpBase[1];
                Point3f p23 = cpBase[2];

                // We already have p12.
                Point3f p11 = Lerp(0.5, p01, p12);
                Point3f p22 = Lerp(0.5, p12, p23);

                // Now elevate to degree 3.
                segCpBezier[0] = p11;
                segCpBezier[1] = Lerp(2.f/3.f, p11, p12);
                segCpBezier[2] = Lerp(1.f/3.f, p12, p22);
                segCpBezier[3] = p22;
            } else {
                // Otherwise we will blossom from p012, p123, p234, and p345
                // to the Bezier control points p222, p223, p233, and p333.
                // https://people.eecs.berkeley.edu/~sequin/CS284/IMGS/cubicbsplinepoints.gif
                Point3f p012 = cpBase[0];
                Point3f p123 = cpBase[1];
                Point3f p234 = cpBase[2];
                Point3f p345 = cpBase[3];

                Point3f p122 = Lerp(2.f/3.f, p012, p123);
                Point3f p223 = Lerp(1.f/3.f, p123, p234);
                Point3f p233 = Lerp(2.f/3.f, p123, p234);
                Point3f p334 = Lerp(1.f/3.f, p234, p345);

                Point3f p222 = Lerp(0.5f, p122, p223);
                Point3f p333 = Lerp(0.5f, p233, p334);

                segCpBezier[0] = p222;
                segCpBezier[1] = p223;
                segCpBezier[2] = p233;
                segCpBezier[3] = p333;
            }
            ++cpBase;
        }

        auto c = CreateCurve(o2w, w2o, reverseOrientation, segCpBezier,
                             Lerp(Float(seg) / Float(nSegments), width0, width1),
                             Lerp(Float(seg + 1) / Float(nSegments), width0, width1),
                             type, n ? &n[seg] : nullptr, sd);
        curves.insert(curves.end(), c.begin(), c.end());
    }
    return curves;
}

}  // namespace pbrt
