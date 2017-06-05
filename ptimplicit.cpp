/*
	DEVELOPER:	MCO

	Implicit path tracing with uniform solid angle sampling and biased path-length termination
*/

#include "pathintegrator.h"

NORI_NAMESPACE_BEGIN

class PtImplicit : public PathIntegrator
{
public:
	PtImplicit(const PropertyList& propList)
	{
		// Get path length
		m_iPathLength = propList.getInteger("pathLenght", 5);	//	Default is 5

		// Get number of samples per pixel
		m_iSampleNo = propList.getInteger("samples", 10);	//	Default is 10

		//	!!!!!!	m_fInvPDF	ATA	!!!!!!
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, int indexX = 0, int indexY = 0) const
	{
		/*
			Follow normal ray tracing procedure until finding an intersection.
			
			* If there is no intersection, return background color
			* If this object is light directly return light color
			* Otherwise, call path tracing integrator, and return calculated result

			Shadow rays, and checking whether the object in shadow or not will be
			handle during the MC process
		*/

		//Check for ray-object intersection
		Intersection its;
		if(!scene->rayIntersect(ray, its))	//	No intersection
		{
			return Color3f(0.f);	// Return background (Black) color
		}

		//Check whether the object is light
		if(its.mesh->isLight())
		{
			// Return color of the light source
			return its.mesh->getBSDF()->getAlbedo();
		}

		//Start Path Tracing
		return tracePath(scene, sampler, its, -ray.d);


	}

private:

	Color3f tracePath(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor(0);

		for(int i = 0; i < m_iSampleNo; i++)
		{
			pixelColor += trace(scene, sampler, firstIntersection, firstRayDirection, 0);
		}

		return pixelColor / m_iSampleNo;
	}



	Color3f trace(const Scene *scene, Sampler *sampler, const Intersection its, const Vector3f firstRayDirection, int pathLength) const
	{
		if( pathLength > m_iPathLength)	// Siniri astigi halde isiga carpamadi
		{
			return Color3f(0);
		}

		//return Color3f(1);

		// 1. Pick a direction from the point
		//Vector3f wi = squareToCosineHemisphere(sampler->next2D()); //squareToUniformSphere(sampler->next2D()); //squareToUniformHemisphere(sampler->next2D()); //squareToCosineHemisphere(sampler->next2D());


		/*
				Formulae, which are used for generating random direction on unit sphere, are taken from
				PBRT 2nd Edition, CH 13 Monte Carlo Integration I: Basic Concepts page 664

				x, y, and z values of W direction is calculated for this task.
				During the calculation random points r1, and r2
				are selected using Sampler.next2D() method of Nori.

				Direction W:

				Direction.z = 1 - 2r1
				Direction.x = sqrt(1 - Direction.z^2) * cos(2 * M_PI * r2)
				Direction.y = sqrt(1 - Direction.z^2) * sin(2 * M_PI * r2)

				PDF:	1 / 4PI

				MC:		1/N	[ E (Li * BRDF(P/M_PI) * n dot W) / PDF ]  ->  4/N	[ E (Li * P * n dot W) ]
		*/

			Point2f dummyP = sampler->next2D();
			
			float r1	=	dummyP.x();
			float r2	=	dummyP.y();

			float z = 1.f - (2.f * r1);

			float sqrtValue		=	sqrtf(std::max(0.f, (1.f - (z * z))));	// r
			float cosSinValue	=	2.f * M_PI * r2;	//phi

			//Calculate direction
			Vector3f wi( sqrtValue * cosf(cosSinValue),	//	X
						 sqrtValue * sinf(cosSinValue),	//	Y
						 z);							//	Z
							  
			wi.normalize();



		// 2. Calculate BRDF * cosQ / pdf at the point P
		
		// Get BRDF
		BSDFQueryRecord bRec(its.toLocal(firstRayDirection));	// Query to get BRDF value of the object
																// For color, and shadow calculations Local coordinates are used

		bRec.wo = its.toLocal(wi);

		Color3f brdf = its.mesh->getBSDF()->eval(bRec);

		// Calculate cosQ

		float theta	= std::max(0.f, its.shFrame.n.dot(its.toWorld(wi)));

		// Get pdf

		float invPDF = (4 * M_PI); //M_PI; //4 * M_PI; //2 * M_PI; //M_PI;

		// Calculate Pi result

		Color3f Pi = brdf * theta * invPDF;

		//Color3f Pi = brdf * theta * invPdf;

		//Color3f Pi = Color3f(1,0,0) * invPdf; //brdf * invPdf;

		// 3. Shout the ray if it hits the light returns Le value of the light
		// else call yourself from the latest point you hit

		// Sample a new ray
		Ray3f sampleRay(its.p, its.toWorld(wi));

		// Shoot the ray, and check the reult for light
		Intersection sampleIts;
		if(scene->rayIntersect(sampleRay, sampleIts))	// It intersects with light
		{
			if(sampleIts.mesh->isLight())
			{
				return Pi * sampleIts.mesh->getBSDF()->getAlbedo();
			}
			else	// No intersection with Light
			{
				//pathLength++;	//Increase path length for next step
				return Pi * trace(scene, sampler, sampleIts, -wi, pathLength + 1);	// Call for next path tracing
			}
		}

		// No intersection
		return Color3f(0.f);

	}

	int m_iPathLength;	// Path length of path tracer
	int m_iSampleNo;	// Number of samples
	//float m_fInvPDF;	// Inverse PDF value
};

NORI_REGISTER_CLASS(PtImplicit, "pathtr");
NORI_NAMESPACE_END