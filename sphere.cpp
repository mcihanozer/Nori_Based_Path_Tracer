/*
Basic sphere object that inharites GeometricShape class.
It contains transform (Position and scale (as redius)),
BRDF (color), and name.
*/

#include "sphere.h"

NORI_NAMESPACE_BEGIN

#pragma optimize("", off)


Sphere::Sphere(const PropertyList &propList)
{
	m_transform	=	propList.getTransform("toWorld", Transform());	// Get sphere transform
	m_name		=	propList.getString("objName", "NoNameForYou");

	m_bIsLight	=	propList.getBoolean("isLight", false);
	m_iId		=	propList.getInteger("id", -1);
	m_fKd		=	propList.getFloat("kd", 0.99);

	m_fRadius = m_transform.getMatrix()(0,0);	// Get radius
	m_p3fCentre = new Point3f(m_transform.getMatrix()(0, 3), m_transform.getMatrix()(1, 3), m_transform.getMatrix()(2, 3));	// Set sphere position

	m_vertexPositions = NULL;

	m_vertexPositions = m_p3fCentre;

	m_vertexCount = 1;
	m_triangleCount = 1;

	m_vertexNormals = NULL; //?????
	m_vertexTexCoords = NULL; //?????
	m_indices = NULL; //?????

}

Sphere::~Sphere()
{
	delete[] m_p3fCentre;
	m_p3fCentre = NULL;

	m_vertexPositions = NULL;
}

//Activate BSDF and re-init Mesh attributes for Sphere class
void Sphere::activate()
{
	if (!m_bsdf)
	{
		// If no material was assigned, instantiate a diffuse BRDF
		m_bsdf = static_cast<BSDF *>(NoriObjectFactory::createInstance("diffuse", PropertyList()));
	}
}

void Sphere::addChild(NoriObject *obj)
{
	switch (obj->getClassType())
	{
	case EBSDF:
		if (m_bsdf)
			throw NoriException("Mesh: tried to register multiple BSDF instances!");
		m_bsdf = static_cast<BSDF *>(obj);
		break;

	default:
		throw NoriException(QString("Mesh::addChild(<%1>) is not supported!").arg(
			classTypeName(obj->getClassType())));
	}
}

//METHODS

// Calculatea point on the sphere
/*
	For Spherical Surface Area calculations formulae, Global Illumination Compendium is used
	(LINK:	http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf	)

	For picking a point on the sphere:

	x = Cx + ( 2r * cos(2*PI*u1) * sqrt( u2(1-u2) )
	y = Cy + ( 2r * sin(2*PI*u1) * sqrt( u2(1-u2) )
	z = Cz + r(1 - 2u2)

	formulae are used, where Cx,y,z is the center of sphere, r is radius of the sphere, and u1,2 are random
	points which are selected using Sampler.next2D() method of Nori.
*/
Point3f Sphere::getPointOn(Point2f randomPoints) const
{
	// Get random values
	float u1 = randomPoints.x();
	float u2 = randomPoints.y();

	// Calculate sqrt, and cos-sin values
	float sqrtVal	=	sqrt(std::max( 0.f, (u2 * (1-u2)) ));
	float cosSinVal	=	2 * M_PI * u1;

	//Calculate and return the point
	return Point3f( m_p3fCentre->x() + ( 2 * m_fRadius * cosf(cosSinVal) * sqrtVal),	// X
					m_p3fCentre->y() + ( 2 * m_fRadius * sinf(cosSinVal) * sqrtVal),	// Y
					m_p3fCentre->z() + ( m_fRadius * (1 - 2 * u2)) );					// Z
}


// Calculates normal value at a given point
/*
	For the calculation, SN = [(xi - xc)/Sr,   (yi - yc)/Sr,   (zi - zc)/Sr]
	formula from Siggraph tutorial is used. Where, SN is the surface normal,
	Xi, Yi, Zi are given point coordinates, Xc, Yc, Zc are coordinates
	of the sphere's center.
*/
Normal3f Sphere::getNormalAt(const Point3f& p) const
{
	float invRadius = 1.f / m_fRadius;

	Normal3f normal( (p.x() - m_p3fCentre->x()) * invRadius,	// X
					 (p.y() - m_p3fCentre->y() ) * invRadius,	// Y
					 (p.z() - m_p3fCentre->z()) * invRadius);	// Z
	normal.normalize();

	return normal;
}

//Calculates Ray-Sphere Intersection
/*
This method uses algebraic solution for calculating ray-sphere intersection.
Algorithm of the method is referenced from Siggraph.

Link: https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
*/
bool Sphere::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const	// Ray-Object intersection
{
	//TODO		Note: If |Rd| = 1 (normalized), then A = 1. So we can compute Sr^2 once.
	//			So with A = 1, the solution of the quadratic equation is
	//	BU HALE GETIR!!!

	//	TODO	EFFICIENT HALE GETIR!!!!!!

	//Get values (For avoiding fucntion call overhead) CONS: Memory consumption
	float Xd = ray.d.x();		float Yd = ray.d.y();		float Zd = ray.d.z();
	float Xo = ray.o.x();		float Yo = ray.o.y();		float Zo = ray.o.z();
	float Xc = m_p3fCentre->x();	float Yc = m_p3fCentre->y();	float Zc = m_p3fCentre->z();	

	//float Xc = m_p3fCentre.x();	float Yc = m_p3fCentre.y();	float Zc = m_p3fCentre.z();	

	//Calculate A value
	// A = Xd^2 + Yd^2 + Zd^2 where Xd, Yd, and Zd are direction of the Ray
	float a =  Xd * Xd + Yd * Yd + Zd * Zd;

	//Calculate B value
	// B = 2 * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc)) where X0, Y0, Z0
	// are origin of the ray, and Xc, Yc, Zc are center of the sphere
	float b = 2 * (Xd * (Xo - Xc) + Yd * (Yo - Yc) + Zd * (Zo - Zc));

	//Calculate C value
	// C = (X0 - Xc)^2 + (Y0 - Yc)^2 + (Z0 - Zc)^2 - Sr^2 where Sr is radius of the sphere
	float c =	(Xo - Xc) * (Xo - Xc) +
		(Yo - Yc) * (Yo - Yc) +
		(Zo - Zc) * (Zo - Zc) - m_fRadius * m_fRadius;

	//Calculate Discriminant
	float disc;
	//if( a == 0 || c == 0 )	Epsilon
	if( a <= Epsilon || c <= Epsilon )
	{
		disc =	b*b;
	}
	else
	{
		disc = (b*b) - (4*a*c);
	}

	//Check whether Discriminant < 0 or not

	//if(disc < 0.f)	// No intersection
	if(disc < Epsilon)
	{

		return false;
	}

	//	std::cout<<"SPHERE RAY INTERSECTION TRUE!\n";

	//Calculate t, and UV values

	//Calculate t using t0, t1 = (- B + (B^2 - 4*C)^1/2) / 2A formula

	//	std::cout<<"DISCRIMINANT!\n";

	float sqrtDisc = std::sqrt(disc);	// Taking square root of the discriminant
	float t0 = ( -b - sqrtDisc ) /	(2*a);

	if(0 <= t0)	//if t0 >= 0, there is no need to calculate t1. Because t0 will be the smallest root
	{
		t = t0;
	}
	else
	{
		t = ( -b + sqrtDisc ) /	(2*a);
	}

	return true;
}

//Get BBox of the Sphere
BoundingBox3f Sphere::getBoundingBox(uint32_t index) const
{
	// Xmin or Xmax			 // Ymin or Ymax			 // Zmin or Zmax
	return BoundingBox3f(Point3f(m_p3fCentre->x() - m_fRadius, m_p3fCentre->y() - m_fRadius, m_p3fCentre->z() - m_fRadius),		//MIN POINTS
		Point3f(m_p3fCentre->x() + m_fRadius, m_p3fCentre->y() + m_fRadius, m_p3fCentre->z() + m_fRadius));	//MAX POINTS
}

QString Sphere::toString() const
{
	return QString(
		"SPHERE[\n"
		"  Transform = %1,\n"
		"  name = \"%1\",\n"
		"  bsdf = %4\n"
		"]")
		.arg(m_transform.toString())
		.arg(m_name)
		.arg(indent(m_bsdf->toString()));
}

NORI_REGISTER_CLASS(Sphere, "sphere");

#pragma optimize("", on)

NORI_NAMESPACE_END