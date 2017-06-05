/*
	DEVELOPER:	MCO

	Implicit path tracing with uniform solid angle sampling and biased path-length termination
*/

#include "pathintegrator.h"

#include "..\sphere.h"	// For being able to cast Light (Mesh) object to Sphere object for colour calculations

NORI_NAMESPACE_BEGIN

class PtExplicit : public PathIntegrator
{
public:
	PtExplicit(const PropertyList& propList)
	{
		m_iSampleNo	=	propList.getInteger("samples", 10);
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

		for(int j = 0; j < m_iSampleNo; j++)
		{
			Vector3f terminator(1);	// For determining path termination
			pixelColor += trace(scene, sampler, firstIntersection, firstRayDirection, terminator);
		}

		return pixelColor / m_iSampleNo;
	}

	Color3f trace(const Scene *scene, Sampler *sampler, const Intersection its, const Vector3f firstRayDirection, Vector3f terminator) const
	{
		// Declarations
		Color3f resultColor = 0;
		const std::vector<Mesh*> lightes = scene->getMeshBasedLights();

		// 1. Calculate Light contrubiton
		for(int i = 0; i < lightes.size(); i++)
		{
			const Sphere* currentLight	=	(Sphere*)lightes[i];
			const int currentLightId	=	currentLight->getId();

			// 1.1 Pick a point on sphere
			Point3f selectedPoint = currentLight->getPointOn(sampler->next2D());

			// 1.2 Generate direction & sample ray

			//Generate direction
			Vector3f direction = selectedPoint - its.p;
			direction.normalize();

			Ray3f sampleRay = Ray3f(its.p, direction);	//Sample ray

			// 1.3 Check for visibility, calculate contrubiton
			// if the light (selected point) is visible
			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts))	//Light, and the selected point is visible
			{
				if(sampleIts.mesh->isLight() && currentLightId == sampleIts.mesh->getId())
				{
					// Calculate contrubiton by "Li * BRDF * cosQ * cosQ' * invPDF * invR2"

					const Color3f	Li			=	currentLight->getBSDF()->getAlbedo();			// Li value for color calculation
					const float		invPDF		=	currentLight->getSurfaceArea();					// inverse PDF of the integrator
					const float		cosTheta	=	std::max(0.f, its.shFrame.n.dot(sampleRay.d));	// Calculate cos(theta): Nintersection dot Wi

					// Calculate BRDF
					BSDFQueryRecord bRec(sampleRay.d, firstRayDirection, EMeasure::ESolidAngle);
					const Color3f brdf = its.mesh->getBSDF()->eval(bRec);

					// Calculae cos(theta') -Light-: Nlight dot - Wi
					Normal3f lightNormal		= currentLight->getNormalAt(selectedPoint);
					const float cosThetaPrime	= std::max(0.f, lightNormal.dot(-sampleRay.d));

					// Calculate inverse distance
					//float		dSquare			= calculateDistanceSquaredBetweenPoints(its.p, selectedPoint) - ( Epsilon * Epsilon);
					//const float inverseDistance = 1.f / dSquare;

					float distanceItoP = calculateDistanceBetweenPoints(its.p, selectedPoint) - Epsilon;
					float dSquare = distanceItoP * distanceItoP;
					float inverseDistance = 1.f / dSquare;

					// Calculate result
					resultColor += Li * brdf * cosTheta * cosThetaPrime * invPDF * inverseDistance;
				}

			} // End of if ray-scene intersection

		} //End of for loop of lightes

		// 2. Check for the path length

		float decider = std::min(0.5f, terminator.y());	// Get variable for path termination checking

		if( unifRand() < decider )	// Don't terminate
		{
			terminator /= decider;	// Change terminator for the next control

			// 2.1 Calculate contribution

			// Pick a direction
			Vector3f wi = squareToCosineHemisphere(sampler->next2D());
			wi.normalize();

			// Calculate BRDF * cosQ / ( Pdf * (1-X) ) for cosine-weighted
			// = BRDF * M_PI / (1-X) 

			// Get BRDF
			BSDFQueryRecord bRec(its.toLocal(firstRayDirection));	// Query to get BRDF value of the object
																	// For color, and shadow calculations Local coordinates are used
			bRec.wo = its.toLocal(wi);
			const Color3f brdf = its.mesh->getBSDF()->eval(bRec);

			// Get invPDF
			const float invPDF = M_PI;

			// Get 1 / (1-X)
			//const float invX = 1 / (1 - decider);
			const float X = (1 - decider);

			// Calculate contribution
			const Color3f colourAtP = brdf * invPDF / X; //* invX;

			// Shout the ray if it hits the light returns Le value of the light
			// else call yourself from the latest point you hit

			// Sample a new ray
			Ray3f sampleRay(its.p, its.toWorld(wi));

			// Shoot the ray, and check the reult for light
			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts))	// It intersects with light
			{
				if(sampleIts.mesh->isLight())
				{
					//resultColor += (colourAtP * sampleIts.mesh->getBSDF()->getAlbedo());
					return resultColor; 
				}
				else	// No intersection with Light
				{
					//pathLength++;	//Increase path length for next step
					resultColor += colourAtP * trace(scene, sampler, sampleIts, -wi, terminator);
				}
			}

		} // End of if no path termination

		return resultColor;
	}

	int m_iSampleNo;	// Number of samples
};

NORI_REGISTER_CLASS(PtExplicit, "pathex");
NORI_NAMESPACE_END