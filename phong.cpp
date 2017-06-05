/*
	DEVELOPER:	MCO

	Modified Phong Reflectance Model implementation

	Using the Modified Phong Reflectance Model for Physically Based Rendering by
	Eric P. Lafortune, and Yves D. Willems is used for MC calculation formulae

	LINK: http://mathinfo.univ-reims.fr/IMG/pdf/Using_the_modified_Phong_reflectance_model_for_Physically_based_rendering_-_Lafortune.pdf )
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <nori\bsdf.h>	// For using BSDFQueryRecord object in the case of intersecting with light source

NORI_NAMESPACE_BEGIN
#pragma optimize("", off)
class Phong : public Integrator
{
public:

	Phong(const PropertyList& propList)
	{
		m_iSampleNo			=	propList.getInteger("sample", 1000);
		//m_fSpecularExponent	=	propList.getFloat("specularEx", 200.f);
	}

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

		// Check for ray-object intersection
		Intersection its;
		if(!scene->rayIntersect(ray, its)){	return	Color3f(0.f);	}	// No intersection: Return back color

		if(its.mesh->isLight()){	return	its.mesh->getBSDF()->getAlbedo();	}	// Intersection with light: Return light color

		return calculatePhong(scene, sampler, its, -ray.d);	// Start Phong Calculation
	}

	QString toString() const
	{
		return QString("PHONG");
	}

private:

	Color3f calculatePhong(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const
	{
		//	Phong Formulate:	diffuse + specular = kd * (1/M_PI) + ks * (n+2 / 2*M_PI) * cons^n(alpha)

		Color3f pixelColor = 0.f;

		// Decide diffuse or specular calculation

		// Get object kd value (Diffuse color)
		//Color3f albedo	= firstIntersection.mesh->getBSDF()->getAlbedo();
		float kd = firstIntersection.mesh->getKd();	// 0.99

		for(int i=0; i < m_iSampleNo; i++)
		{
			float nextStep = unifRand();	// Choose u value for deciding next calculation

			if(nextStep < kd)	// Diffuse calculations
			{
				/*
					Phong Diffuse Formulate:	kd * (1/M_PI)

					PDF:	cosQ / M_PI

					Direction wi:	wi(Q, Phi) = ( arccos(sqrt(u1)), 2*M_PI*u2 )

					MC:		Li * kd
				*/


				// 1. Get direction
				Vector3f direction = getCosineWeightedDirection(sampler->next2D());

				// 2. Sample ray, and check visibility
				Ray3f sampleRay(firstIntersection.p, firstIntersection.toWorld(direction));

				Intersection sampleIts;
				if(scene->rayIntersect(sampleRay, sampleIts))	// If it hits an object
				{
					if(sampleIts.mesh->isLight())	// If the object is light
					{
						// 3. If the point is visible, calculate contribution
						Color3f Li = sampleIts.mesh->getBSDF()->getAlbedo();

						pixelColor +=	Li * firstIntersection.mesh->getBSDF()->getAlbedo(); //albedo;
					}
				}
			}
			else if(nextStep < 1)	// Specular calculations
			//else
			{
				/*
					Phong Specular Formulate:	ks * (n+2 / 2*M_PI) * cons^n(alpha)

					PDF:	(n + 1 / 2*M_PI) * cons^n(alpha)

					Direction wi:	wi(Q, Phi) = ( arccos( u1^(1/n + 1) ), 2*M_PI*u2 )

					MC:		Li * ks * (n + 2 / n + 1) * cosQ
														cosQ = N * wi
				*/

				// 1. Get direction
				float specularExponent = firstIntersection.mesh->getSpecularEx();

				Point2f dummyP = sampler->next2D();

				float r1 = dummyP.x();
				float r2 = dummyP.y();

				float specularExValue = (1 / (specularExponent + 1));
				float theta	=	acosf( pow(r1, specularExValue) );
				float phi	=	2 * M_PI * r2;

				Vector3f direction = Vector3f( sinf(theta) * cosf(phi),		// X
											   sinf(theta) * sinf(phi),		// Y
											   cosf(theta));				// Z
				direction.normalize();

				// 2. Sample ray, and check visibility
				Ray3f sampleRay(firstIntersection.p, firstIntersection.toWorld(direction));

				Intersection sampleIts;
				if(scene->rayIntersect(sampleRay, sampleIts))	// If it hits an object
				{
					if(sampleIts.mesh->isLight())	// If the object is light
					{
						// 3. If the point is visible, calculate contribution

						Vector3f reflectionRay = reflect(firstRayDirection, firstIntersection.shFrame.n);
						float alpha = sampleRay.d.dot( reflectionRay);

						if(alpha > 0 )
						{
							Color3f Li = sampleIts.mesh->getBSDF()->getAlbedo();	// Get Li;

							specularExValue += 1;	// Get (n + 2 / n + 1) = 1 + (1 / n + 1)

							float theta = std::max(0.f, firstIntersection.shFrame.n.dot(sampleRay.d));

							pixelColor += Li * firstIntersection.mesh->getBSDF()->getAlbedo() * specularExValue * theta;
						}
						
					}

				}

			}
			else{	return Color3f(0.f);	}	// No contribution
		}

		return pixelColor / m_iSampleNo;
	}

	int		m_iSampleNo;
	//float	m_fSpecularExponent;
};
#pragma optimize("", on)
NORI_REGISTER_CLASS(Phong, "phong");
NORI_NAMESPACE_END