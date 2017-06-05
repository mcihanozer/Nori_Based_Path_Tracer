/*
	DEVELOPER:	MCO

	Diffuse Monte Carlo Integration
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <nori\bsdf.h>	// For using BSDFQueryRecord object in the case of intersecting with light source
#include "..\sphere.h"	// For being able to cast Light (Mesh) object to Sphere object for explicit MC calculations

NORI_NAMESPACE_BEGIN
#pragma optimize("", off)
class MonteCarlo : public Integrator
{
public:
	//CONSTRUCTOR

	//	Initialize MC Integrator type

	//	----	IMPLICIT	---	///
	//	0: Spherical
	//	1: Hemispherical
	//	2: Cosine-weighted

	//	----	EXPLICIT	---	///
	//	3: Surface Area
	//	4: Importance sampling according to solid angle

	MonteCarlo(const PropertyList& propList)
	{
		// Get integrator type
		m_eIntegratorType = propList.getInteger("type", HEMISPHERICAL_MC);	//	Default is Hemispherical

		if(m_eIntegratorType < SPHERICAL_MC || m_eIntegratorType > SOLIDANGLE_MC)	// If the value assigned is not spherical, hemispherical etc.,
																						// it is assigned hemispherical
		{
			m_eIntegratorType = HEMISPHERICAL_MC;
		}

		// Get number of samples per pixel
		m_iSampleNo = propList.getInteger("samples", 10);	//	Default is 10
	}

	//ABSTRACT CLASS DECLARATIONS

	//	Calculates the color of the pixel using N samples, and spherical, hemispherical or cosine-weighted approaches
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, int indexX, int indexY ) const
	{
		/*
			Follow normal ray tracing procedure until finding an intersection.
			
			* If there is no intersection, return background color
			* If this object is light directly return light color
			* Otherwise, call MC integration, and return calculated result

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
			//return Color3f(1,0,0);

			return its.mesh->getBSDF()->getAlbedo();//->eval(bRec);
		}

		//Start MC intergration
		return CalculateMC(scene, sampler, its, -ray.d);
	}

	QString toString() const
	{
		return QString("Monte Carlo Integration[TYPE=%1]").arg(m_eIntegratorType);
	}

private:
	//	Call suitable MC integration method, and performs Diffuse MC Integration
	Color3f CalculateMC(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor = 0;

		switch(m_eIntegratorType)
		{
		case SPHERICAL_MC:

			//call spherical calculation
			return calculateSphericalMC(scene, sampler, firstIntersection, firstRayDirection);

		case HEMISPHERICAL_MC:

			//call hemispherical calculation
			return calculateHemisphericalMC(scene, sampler, firstIntersection, firstRayDirection);

		case COSINEWEIGHTED_MC:

			//call cosine-weighted calculation
			return calculateCosineMC(scene, sampler, firstIntersection, firstRayDirection);

		case SURFACEAREA_MC:

			//Call surface area calculation (Explicit)
			return calculateSurfaceAreaMC(scene, sampler, firstIntersection, firstRayDirection);

		case SOLIDANGLE_MC:

			// Call solid angle calculation (Explicit)
			return calculateSolidAngleMC(scene, sampler, firstIntersection, firstRayDirection);

		default:
			std::cout<<"\n\n ERROR at MC SWITCH!!!\n\n";
		}

		return pixelColor;
	}

	//Calculates implicit spherical MC Integration
	Color3f calculateSphericalMC(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor = 0;
		float invPDF = (4 * M_PI);

		for(int i=0; i < m_iSampleNo; i++)
		{
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
			Vector3f direction( sqrtValue * cosf(cosSinValue),	//	X
								sqrtValue * sinf(cosSinValue),	//	Y
								z								//	Z
							  );
			direction.normalize();

			//Create sample ray
			Ray3f sampleRay(firstIntersection.p, firstIntersection.toWorld(direction));	//Converts the ray direction to World/Global coordinates
																						// For checking intersection with light source, World/Global coordinates are used

			//Fire the ray, and check the result

			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts))	// If it hits an object
			{
				if(sampleIts.mesh->isLight())	// And the object is light
				{
					BSDFQueryRecord bRec(firstIntersection.toLocal(firstRayDirection));	// Query to get BRDF value of the object
																						// For color, and shadow calculations Local coordinates are used
					bRec.wo = firstIntersection.toLocal(direction);

					// Calculate color of the sample

					Color3f Li		= sampleIts.mesh->getBSDF()->getAlbedo();
					Color3f brdf	= firstIntersection.mesh->getBSDF()->eval(bRec);
					float theta		= std::max(0.f, firstIntersection.shFrame.n.dot(sampleRay.d));

					pixelColor		+= Li * brdf * theta * invPDF;
				}
			}

		} //End of for

		return pixelColor / m_iSampleNo;
	}

	//Calculates implicit hemispherical MC Integration
	Color3f calculateHemisphericalMC(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor = 0;
		float invPDF = (2 * M_PI);	// PDF:	1 / 2PI calculation

		for(int i=0; i < m_iSampleNo; i++)
		{
			/*
				Formulae, which are used for generating random direction on unit hemisphere proportional
				to solid angle, are taken from Global Illumination Compendium
				(LINK:	http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf	)
				 
				x, y, and z values of W direction is calculated for this task.
				During the calculation random points r1, and r2
				are selected using Sampler.next2D() method of Nori.

				Direction W:

				Direction.x = cos(2*M_PI*r1) * sqrt(1 - r2^2)
				Direction.y = sin(2*M_PI*r1) * sqrt(1 - r2^2)
				Direction.z = r2

				PDF:	1 / 2PI

				MC:		1/N	[ E (Li * BRDF(P/M_PI) * n dot W) / PDF ]  ->  2/N	[ E (Li * P * n dot W) ]
			*/

			// Pick a random hemispherical direction
			Point2f dummyP = sampler->next2D();
			
			float r1	=	dummyP.x();
			float r2	=	dummyP.y();

			float sqrtValue		=	sqrtf(std::max(0.f, (1 - (r2*r2))));
			float cosSinValue	=	2 * M_PI * r1;

			//Calculate direction
			Vector3f direction( cos(cosSinValue) * sqrtValue,	//	X
								sin(cosSinValue) * sqrtValue,	//	Y
								r2								//	Z
							  );
			direction.normalize();

			//Create sample ray
			Ray3f sampleRay(firstIntersection.p, firstIntersection.toWorld(direction)); //Converts the ray direction to World/Global coordinates
																						// For checking intersection with light source, World/Global coordinates are used
			//Fire the ray, and check the result
			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts))	// If it hits an object
			{
				if(sampleIts.mesh->isLight())	// If the object is light
				{
					BSDFQueryRecord bRec(firstIntersection.toLocal(firstRayDirection));	// Query to get BRDF value of the object
																						// For color, and shadow calculations Local coordinates are used
					bRec.wo = firstIntersection.toLocal(direction);

					//Calculate the color of the sample

					Color3f Li		=	sampleIts.mesh->getBSDF()->getAlbedo();
					Color3f brdf	=	firstIntersection.mesh->getBSDF()->eval(bRec);
					float theta		=	std::max(0.f, firstIntersection.shFrame.n.dot(sampleRay.d));

					pixelColor		+=	Li * brdf * theta * invPDF;
				}
			}

		} // End of for

		return pixelColor / m_iSampleNo;
	}

	//Calculates implicit Cosine-weighted MC Integration
	Color3f calculateCosineMC(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor = 0;
		float invPDF = M_PI;	// PDF:	cos(theta) / PI but, cos(theta) value is eleminated

		for(int i=0; i < m_iSampleNo; i++)
		{
			/*
				Formulae of IFT 6095 slides are used to implement for picking direction respect to
				unit hemisphere proportional to cosine-weighted solid angle.

				For direction picking, formula on Techniques: illumination globale slide (page 48),
				and for converting to Cartesian coordinates, formula on Th¨¦orie: illumination directe
				slide (page 39) are used.

				...

				PDF:	cos(theta) / PI

				MC:		1/N	[ E (Li * BRDF(P/M_PI) * n dot W) / PDF ]  ->  1/N	[ E (Li * P) ]
			*/

			Point2f dummyP = sampler->next2D();

			float r1 = dummyP.x();
			float r2 = dummyP.y();

			float theta = acosf( sqrt(r1) );
			float phi	= 2 * M_PI * r2;

			Vector3f direction = Vector3f( sinf(theta) * cosf(phi),		// X
										   sinf(theta) * sinf(phi),		// Y
										   cosf(theta));				// Z


			//Vector3f direction = squareToCosineHemisphere(sampler->next2D());
			direction.normalize();

			//Create sample ray
			Ray3f sampleRay(firstIntersection.p, firstIntersection.toWorld(direction));

			//Fire the ray and check the result
			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts)) // If it hits an object
			{
				if(sampleIts.mesh->isLight()) // If the object is light
				{
					BSDFQueryRecord bRec(firstIntersection.toLocal(firstRayDirection)); // Query to get BRDF value of the object
																						// For color, and shadow calculations Local coordinates are used

					bRec.wo = firstIntersection.toLocal(direction);

					//Calculate the color of the sample

					Color3f Li		= sampleIts.mesh->getBSDF()->getAlbedo();
					Color3f brdf	= firstIntersection.mesh->getBSDF()->eval(bRec);
					
					pixelColor		+= Li * brdf * invPDF;
				}
			}

		} // End of for

		return pixelColor / m_iSampleNo;
	}

	Color3f calculateSurfaceAreaMC(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor = 0;
		const std::vector<Mesh*> lightMeshes = scene->getMeshBasedLights();

		//For all light sources, sample the light and calculate the result for the point
		for(int i = 0; i < lightMeshes.size(); i++)
		{
			pixelColor += sampleSurfaceArea(scene, sampler, firstIntersection, firstRayDirection, lightMeshes[i]);
		}

		return pixelColor;

	}//End of calculateSurfaceAreaMC()

	Color3f sampleSurfaceArea(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection, const Mesh* light) const
	{
		/*
				For Spherical Surface Area calculations formulae, Global Illumination Compendium is used
				(LINK:	http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf	)

				For picking a point on the sphere:

				x = Cx + ( 2r * cos(2*PI*u1) * sqrt( u2(1-u2) )
				y = Cy + ( 2r * sin(2*PI*u1) * sqrt( u2(1-u2) )
				z = Cz + r(1 - 2u2)

				formulae are used, where Cx,y,z is the center of sphere, r is radius of the sphere, and u1,2 are random
				points which are selected using Sampler.next2D() method of Nori.

				PDF : 1 / 4*PI*r^2 where r is radius of the sphere (1 / Surface area of sphere)

				MC:		1/N	[ E (Li * BRDF(P/M_PI) * n dot W * nL dot WL) / PDF * d^2 ]
				
				, where nL dor wL is the cosine value between the selected direction, and light source normal, and d^2 is
				the distance between intersection point, and selected point on the light


		*/

		//Init variables
		Color3f lightContribution = 0;

		const Sphere*	currentLight	=	(Sphere*)light;
		const Point3f	center			=	currentLight->getCenter();				// Center for sampling point
		const float		radius			=	currentLight->getRadius();				// Radius for sampling point
		const float		invPDF			=	currentLight->getSurfaceArea();			// inverse PDF of the integrator
		const Color3f	Li				=	currentLight->getBSDF()->getAlbedo();	// Li value for color calculation
		const int		currentLightId	=	currentLight->getId();					// ID of the light for checking visibility

		//Sample the light
		for(int i = 0; i < m_iSampleNo; i++)
		{
			// 1. Pick point

			// Selecting random values
			Point2f dummyP = sampler->next2D();
			float u1 = dummyP.x();
			float u2 = dummyP.y();

			// Calculate sqrt, and cos-sin values
			float sqrtVal	=	sqrt(std::max( 0.f, (u2 * (1-u2)) ));
			float cosSinVal	=	2 * M_PI * u1;

			//Calculate point
			Point3f selectedPoint = Point3f( center.x() + ( 2 * radius * cosf(cosSinVal) * sqrtVal),	// X
											 center.y() + ( 2 * radius * sinf(cosSinVal) * sqrtVal),	// Y
											 center.z() + ( radius * (1 - 2 * u2)) );					// Z

			// 2. Calculate direction (world coordinate)
			Vector3f direction = selectedPoint - firstIntersection.p;
			direction.normalize();

			// 3. Generate sample ray
			Ray3f sampleRay = Ray3f(firstIntersection.p, direction); //, Epsilon, distanceItoP);

			// 4. Check whether the point selected on the light is visible for the intersection point
			// If so, calculate contribution
			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts)) //Light, and the selected point is visible
			{
				if(sampleIts.mesh->isLight() && sampleIts.mesh->getId() == currentLightId)
				{
					// Calculate MC:		Li * BRDF * cosQ * cosQ' * invPDF * invR2

					// Calculate BRDF
					BSDFQueryRecord bRec(sampleRay.d, -firstRayDirection, EMeasure::ESolidAngle);
					Color3f brdf = firstIntersection.mesh->getBSDF()->eval(bRec);

					// Calculate cos(theta): Nintersection dot Wi
					float cosTheta = std::max(0.f, firstIntersection.shFrame.n.dot(sampleRay.d));

					// Calculae cos(theta') -Light-: Nlight dot - Wi
					Normal3f lightNormal = currentLight->getNormalAt(selectedPoint);
					float cosThetaPrime = std::max(0.f, lightNormal.dot(-sampleRay.d));

					// Calculate inverse distance
					float distanceItoP = calculateDistanceBetweenPoints(firstIntersection.p, selectedPoint) - Epsilon;
					float dSquare = distanceItoP * distanceItoP;
					float inverseDistance = 1.f / dSquare;

					// Calculate result
					lightContribution += Li * brdf * cosTheta * cosThetaPrime * invPDF * inverseDistance;
				}
			}

		}

		return lightContribution /= m_iSampleNo;
	}

	Color3f calculateSolidAngleMC(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		Color3f pixelColor = 0;
		const std::vector<Mesh*> lightMeshes = scene->getMeshBasedLights();

		// For all light sources, sample the light, and calculate the result for the point
		for(int i = 0; i < lightMeshes.size(); i++)
		{
			pixelColor += sampleSolidAngle(scene, sampler, firstIntersection, firstRayDirection, lightMeshes[i]);
		}

		return pixelColor;
	}

	Color3f sampleSolidAngle(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection, const Mesh* light) const
	{
		/*
			Spherical light importance sampling implementation according to solid angle.
			
			For coding this function, PBRT 2nd edition book, Section 14.6 Sampling Light Sources
			[720,722] pages, and PBRT v2 source code are referenced.
		*/

		// Init variables
		Color3f lightContribution;

		const Sphere*	currentLight	=	(Sphere*)light;
		const Point3f	center			=	currentLight->getCenter();				// Center for sampling point
		const float		radius			=	currentLight->getRadius();				// Radius for sampling point
		const Color3f	Li				=	currentLight->getBSDF()->getAlbedo();	// Li value for color calculation
		const int		currentLightId	=	currentLight->getId();					// ID of the light for checking visibility

		// Sample Light
		for(int i = 0; i < m_iSampleNo; i++)
		{
			// 1. Get direction: Needs random values, center of the sphere (light), normal at the intersection point,
			// cosQmax (needs radius of the sphere (light), and the coordinate system (needs vector between intersection
			// point, and center of the sphere (light) for aligning the cone with the light

			// Get random values
			Point2f dummyP = sampler->next2D();
			float u1 = dummyP.x();
			float u2 = dummyP.y();

			// Calculate variables for coordinate system calculation
			Vector3f vecInttoCen = center - firstIntersection.p;	// Vector between intersection point, and center of the sphere (light)
			vecInttoCen.normalize();
			Vector3f wcX, wcY;
			coordinateSystem(vecInttoCen, wcX, wcY);	// Create a new coordinate system

			// Calculate cosQmax
			float sinTmax2		=	radius * radius / calculateDistanceSquaredBetweenPoints(center, firstIntersection.p);
			float cosThetaMax	=	sqrtf(std::max(0.f, 1.f -  sinTmax2));

			//Get direction
			Vector3f rayDirection = uniformSampleCone(u1, u2, cosThetaMax, wcX, wcY, vecInttoCen);
			rayDirection.normalize();

			// 2. Sample the ray
			Ray3f sampleRay(firstIntersection.p, rayDirection);

			// 3. Calculate visibility
			Intersection sampleIts;
			if(scene->rayIntersect(sampleRay, sampleIts))
			{
				if(sampleIts.mesh->isLight() && sampleIts.mesh->getId() == currentLightId)
				{

					// 4. Calculate contrution Li * brdf * cosQ / pdf

					// Calculate BRDF
					BSDFQueryRecord bRec(sampleRay.d, -firstRayDirection, EMeasure::ESolidAngle);
					Color3f brdf = firstIntersection.mesh->getBSDF()->eval(bRec);

					// Calculate cos(theta): Nintersection dot Wi
					float cosTheta = std::max(0.f, firstIntersection.shFrame.n.dot(sampleRay.d));

					float PDF = uniformConePdf(cosThetaMax);

					lightContribution += ( (Li * brdf * cosTheta) / PDF );
				}
			}
		}

		return lightContribution /= m_iSampleNo;
	}

	int m_eIntegratorType;	// Keeps MC integrator type (Spherical, hemispherical etc.)
	int m_iSampleNo;
};
#pragma optimize("", on)
NORI_REGISTER_CLASS(MonteCarlo, "montecarlo");
NORI_NAMESPACE_END