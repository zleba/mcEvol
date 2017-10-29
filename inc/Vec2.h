
class Vec2;

class Vec2 {

	public:

		Vec2(Double xIn = 0., Double yIn = 0.)
			: px(xIn), py(yIn) { }

		Vec2(const Vec2& v) : px(v.px), py(v.py) { }

		void setRPhi(Double  r, Double phi )
		{ px=r*cos(phi); py=r*sin(phi); }
		void setXY(Double  xx, Double yy )
		{ px=xx; py=yy; }



		void x(Double xIn) {px = xIn;}
		void y(Double yIn) {py = yIn;}

		Double x() const {return px;}
		Double y() const {return py;}
		Double pX() const {return px;}
		Double pY() const {return py;}

		Double norm()  const {return sqrt(px*px+py*py);}
		Double norm2() const {return (px*px+py*py);}

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
		friend Double cross(const Vec2& v1, const Vec2& v2);


	private:
		Double px, py;


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


inline Double cross(const Vec2& v1, const Vec2& v2)
{return  v1.px*v2.py - v1.py*v2.px;}

