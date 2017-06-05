/*
	Base class for Light source.
*/

#if !defined(__LIGHTSOURCE_H)
#define __LIGHTSOURCE_H

#include <nori/object.h>
#include<nori\color.h>

NORI_NAMESPACE_BEGIN
/*
	DEVELOPER:	MCO
	
	Base class for Light source.
	
	It might be inherited from Mesh class in the future because,
	for some light sources, it might be in the KD-tree.
*/

class LightSource : public NoriObject
{
public:
	//CONSTRUCTORS
	LightSource()	:	m_c3fLightIntensity(0.f)	{}

	LightSource(float c)	:	m_c3fLightIntensity(c)	{}

	LightSource(float r, float g, float b)	:	m_c3fLightIntensity(r, g, b)	{}

	//ABSTRACT CLASS DECLARATIONS
	virtual EClassType getClassType()	const {	return ELight; }

	virtual Color3f getIntensity()		const {	return m_c3fLightIntensity; }

	//virtual const Point3f getPosition()	{ return Point3f(0);	}

	// Light direction calcualtion method. Because of being different for different
	// light sources, it's an virtual abstract class. In this way, we can eliminate
	// if-else case for every light source.
	//
	// it takes intersection point as parameter which is needed only for point light,
	// area light(?) etc.
	virtual const Vector3f getDirection(Point3f intersectionPosition) = 0;

	//ABSTRACT METHOD DEFINITION
	/*
		Takes surface normal and calculates cosine value between light direction
		and surface value.

		The calculation (except dot product part) is different for directional and
		point light sources
	*/
	//virtual float getCosValue(const Normal3f surfaceNormal) = 0;

protected:
	Color3f m_c3fLightIntensity;
};

NORI_NAMESPACE_END

#endif