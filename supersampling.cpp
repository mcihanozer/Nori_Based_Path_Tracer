/*
	Developer:	MCO

	This file contains SuperSampling and SuperSamplingRandom classes.
	
	SuperSampling is uniform implementation of Super Sampling anti-aliasing
	approach. Currently, it always shoots 4 rays to one pixel in	*		*
																		*
																	*		*
	orientation.

	SuperSamplingRandom is stochastic implementation of Super Sampling anti-aliasing
	approach. It can shoot any number of rays in a random fashion. It uses std::rand()
	function of C++ for selecting randim positions.
*/

#include <nori/object.h>

#include "..\antialiasing.h"

NORI_NAMESPACE_BEGIN

class SuperSampling : public AntiAliasing
{
public:
	SuperSampling(const PropertyList &propList)
	{
		isRandom	= false;
		aaCount		= 4; //arrangeCount(propList.getInteger("sampleCount", 4));

		aaPoints.push_back(Point3f(0.1666666f, 0.1666666f, 0));	//Position #1
		aaPoints.push_back(Point3f(0.8333333f, 0.1666666f, 0));	//Position #2
		aaPoints.push_back(Point3f(0.1666666f, 0.8333333f, 0));	//Position #3
		aaPoints.push_back(Point3f(0.8333333f, 0.8333333f, 0));	//Position #4
	}

	virtual QString toString() const
	{
		return QString(
			"Anti-Aliasing[\n"
			"  Count = %3,\n"
			"  Name = \n"
			"]")
		.arg(aaCount)
		.arg("Uniform");
	}
};
NORI_REGISTER_CLASS(SuperSampling, "uniform");

class SuperSamplingRandom : public AntiAliasing
{
public:
	/*SuperSamplingRandon(int count)
	{
		isRandom = true;
		aaCount = count;
	}*/

	SuperSamplingRandom(const PropertyList &propList)
	{
		isRandom	= true;
		aaCount		= propList.getInteger("sampleCount", 4); //arrangeCount(propList.getInteger("sampleCount", 4));
	}

	virtual QString toString() const
	{
		return QString(
			"Anti-Aliasing[\n"
			"  Count = %3,\n"
			"  Name = \n"
			"]")
		.arg(aaCount)
		.arg("Stochastic");
	}
};
NORI_REGISTER_CLASS(SuperSamplingRandom, "stochastic");
NORI_NAMESPACE_END