/*
	DEVELOPER:	MCO

	Directional Light source. It inherits LightSource class and
	adds light direction to the implementation.
*/


#include <nori/vector.h>

#include "lightsource.h"

NORI_NAMESPACE_BEGIN
#pragma optimize("", off)
class DirectionalLight	: public LightSource
{
public:
	//CONSTRUCTORS
	DirectionalLight() :	LightSource(0)
	{
		m_v3fDirection.x() = 0.f;
		m_v3fDirection.y() = 0.f;
		m_v3fDirection.z() = 0.f;
	}

	DirectionalLight(const PropertyList &propList)
	{
		m_v3fDirection		= propList.getVector("direction");	//Get Direction vector

		m_v3fDirection.normalize();

		Vector3f dumy = propList.getVector("intensity");

		m_c3fLightIntensity = Color3f( dumy.x(), dumy.y(), dumy.z());  //propList.getVector("intensity");	// Get Color
	}

	//ABSTRACT CLASS DECLARATIONS
	/*
	Calculates angle between intersected surface normal and light direction
	and returns the value.

	For the calculation, 'N dot L' formula is used. there, N stads for surface
	normal, L stands for direction of the light.
	*/
	//virtual const Vector3f getDirection() { return m_v3fDirection; }

	virtual const Vector3f getDirection(Point3f intersectionPosition) { return m_v3fDirection; }

	/*virtual float getCosValue(const Normal3f surfaceNormal)
	{
		return surfaceNormal.dot(m_v3fDirection);
	}*/

	QString toString() const
	{
		return QString(
			"Light Soruce[\n"
			"  name = \"%1\",\n"
			"  Direction = %2,\n"
			"  Intensift = %3,\n"
			"]")
			.arg("Directional Light")
			.arg(m_v3fDirection.toString())
			.arg(m_c3fLightIntensity.toString());
	}


private:
	Vector3f m_v3fDirection;
};
#pragma optimize("", on)
NORI_REGISTER_CLASS(DirectionalLight, "directional");
NORI_NAMESPACE_END