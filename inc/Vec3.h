
class Vec4;

class Vec4 {

	public:

		Vec4(Double eIn=0, Double xIn = 0., Double yIn = 0., Double zIn=0)
			: e(eIn), px(xIn), py(yIn), pz(zIn) {}

		Vec4(const Vec4& v) : e(v.e), px(v.px), py(v.py), pz(v.pz) { }

		void setRPhi(Double  r, Double phi)
		{ px=r*cos(phi); py=r*sin(phi); }
		void setXY(Double  xx, Double yy )
		{ px=xx; py=yy; }


        void AddPt(double z, double ) {
            Vec4 pNew = z*this;
            //(1,0,0)
            Vec4 pNew1(0, -pz, py);
            Vec4 pNew2(-pz, 0, px);


        }

		void x(Double xIn) {px = xIn;}
		void y(Double yIn) {py = yIn;}

		Double x() const {return px;}
		Double y() const {return py;}
		Double pX() const {return px;}
		Double pY() const {return py;}

		Double norm()  const {return sqrt(px*px+py*py+pz*pz);}
		Double mass2() const {return (e*e-px*px-py*py-pz*pz);}

		Vec2 operator-()
		{ return Vec2(-px,-py);}
		Vec2& operator+=(const Vec2& v)
		{px += v.px; py += v.py; return *this;}
		Vec2& operator-=(const Vec2& v)
		{px -= v.px; py -= v.py; return *this;}
		Vec2& operator*=(Double f)
		{px *= f; py *= f; return *this;}
		Vec2& operator/=(Double f)
		{px /= f; py /= f; return *this;}

		// Operator overloading with friends
		friend Vec2 operator+(const Vec2& v1, const Vec2& v2);
		friend Vec2 operator-(const Vec2& v1, const Vec2& v2);
		friend Vec2 operator*(Double f, const Vec2& v1);
		friend Vec2 operator*(const Vec2& v1, Double f);
		friend Vec2 operator/(const Vec2& v1, Double f);
		friend Double operator*(const Vec2& v1, const Vec2& v2);
		friend Double cross(const Vec4& v1, const Vec4& v2);


	private:
		Double e, px, py, pz;


};

inline Vec2 operator+(const Vec2& v1, const Vec2& v2)
{Vec2 v = v1 ; return v += v2;}

inline Vec2 operator-(const Vec2& v1, const Vec2& v2)
{Vec2 v = v1 ; return v -= v2;}

inline Vec2 operator*(Double f, const Vec2& v1)
{Vec2 v = v1; return v *= f;}

inline Vec2 operator*(const Vec2& v1, Double f)
{Vec2 v = v1; return v *= f;}

inline Vec2 operator/(const Vec2& v1, Double f)
{Vec2 v = v1; return v /= f;}

inline Double operator*(const Vec2& v1, const Vec2& v2)
{return  v1.px*v2.px + v1.py*v2.py;}


inline vec4 cross(const Vec4& v1, const Vec4& v2)
{
    return vec4(0, v1.py*v2.pz - v1.pz*v2py, v1.pz*v2.px - v1.px*v2.pz, v1.px*v2.py - v1.py*v2.px);
}

