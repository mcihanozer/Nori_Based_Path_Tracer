/*
	Developer:	MCO

	Base class for anti-aliasing. It is used for creating an high
	level approaches. Camera object (perspective.cpp) is used
	AntiAliasing object for storing the anti-aliasing object and
	using it in diffuse illumunation in do.cpp file.

	Currently, only Super Sampling approach is supported. Blue noise
	or jittering algorithm will be implemented in the future.
*/

#if !defined(__ANTIALIASING_H)
#define __ANTIALIASING_H

#include<vector>

#include<nori\vector.h>

NORI_NAMESPACE_BEGIN

class AntiAliasing :  public NoriObject
{
public:
	virtual EClassType getClassType()	const	{	return EAAliasing;	}
	virtual QString toString()			const	{	return QString("Base Anti-aliasing"); }

	bool isRandom;
	int aaCount;

	std::vector<Point3f> aaPoints;	// Keeps number of anti-aliasing rays

protected:
	int arrangeCount(int count)
	{
		if(count <= 8)
		{
			return 4;
		}
		else if(count <= 32)
		{
			return 16;
		}
		else
		{
			return 64;
		}
	}
};

NORI_NAMESPACE_END

#endif