/*
	Developer:	MCO

	Basic sphere object that inharites GeometricShape class.
	It contains transform (Position and scale (as redius)),
	BRDF (color), and name.
*/

#include <nori/bbox.h>

#include "geometricobjects.h"

NORI_NAMESPACE_BEGIN

class Sphere : public GeometricShape
{
public:

	Sphere(const PropertyList &propList);
	virtual ~Sphere();

	void activate();
	void addChild(NoriObject *obj);


	//GETTER
	inline const float		getRadius()			const	{	return m_fRadius;	 }
	inline const Point3f	getCenter()			const	{	return *m_p3fCentre; }
	inline const float		getSurfaceArea()	const	{	return 4*M_PI*m_fRadius*m_fRadius;	}

	Point3f getPointOn(Point2f randomPoints) const;

	Normal3f getNormalAt(const Point3f& p) const;

	bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const;

	BoundingBox3f getBoundingBox(uint32_t index) const;

	BoundingBox3f getClippedBoundingBox(uint32_t index, const BoundingBox3f &bbox) const
	{
		return getBoundingBox(index);
	}

	QString toString() const;

private:
	
	float		m_fRadius;
	Point3f		*m_p3fCentre;

};

NORI_NAMESPACE_END