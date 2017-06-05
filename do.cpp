/*
	DEVELOPER:	MCO

	Diffuse rendering class
*/

#include <time.h> 

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <nori\camera.h>

#include "..\geometricobjects.h"
#include "..\antialiasing.h"

NORI_NAMESPACE_BEGIN

 /*
	Diffuse Shading class for Nori. It inhertis Integrator Class and
	implements Li() method which calculates the color of the pixel.

	For the calculation, Diffuse shading formula
	I = IaKa + IpKd(N dot L) is used.
	
	There, I stands for intenstiy (color) of the pixel, Ia for Ambient
	intensity and Ka for Ambient coefficient. IaKa together defines
	Ambient color.
	
	Ip stands for intensity of the light source and Kd for diffuse
	reflection coeffient of the object that means BRDF of the object
	and the value is P/PI for Diffuse Shading.
	
	(N dot L) defines dot product between surface normal (N) and light
	direction vector (L). The result is the cosine value (angle) between
	these two vectors.
	
	IpKd(N dot L) stands for Diffuse color.
*/
#pragma optimize("", off)
class DiffuseOcclusion : public Integrator
{
public:
	//CONSTRUCTORS
	DiffuseOcclusion() : m_fLength(0.f)	{}

	DiffuseOcclusion(const PropertyList& propList)
	{
		m_fLength = propList.getFloat("length", 1e4f);	// Ray length of the ambient occlusion queries expressed relative to the scene size
	}

	//ABSTRACT CLASS DECLARATIONS
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, int indexX, int indexY ) const
	{
		///		!!!		BLACK = 0	WHITE = 1		!!!		///

		Color3f pixelColor = 0.f;
		const AntiAliasing* antialiaser = scene->getCamera()->getAntialiaser();
		const Transform cameraTransform = scene->getCamera()->getCameraTransform();

		//Send anti-aliasing rays and calculate colors
		for(int i = 0; i < antialiaser->aaCount; i++)
		{
			//Generate anti-aliasing rays

			//TODO	Her iterationda 'isRandom?' kontrolu yapmak gereksiz. Buna daha uygun bir yontem bul.
			//		Rayleri fordan once generate edip sonra for icinde gondermek bir cozum olabilir!

			//Create next ray's point on the image plane for ray direction calculations
			Point2f pixelSample;
			if( antialiaser->isRandom)
			{
				pixelSample = Point2f((std::rand() / (double)RAND_MAX) + indexX, (std::rand() / (double)RAND_MAX) + indexY);	//	stochastic
			}
			else
			{
				pixelSample = Point2f(antialiaser->aaPoints[i].x() + indexX, antialiaser->aaPoints[i].y() + indexY);	//	uniform
			}
			
			Point2f apertureSample = sampler->next2D();
			//Sample an anti-aliasing ray from the camera
			Ray3f aaRay;
			Color3f value = scene->getCamera()->sampleRay(aaRay, pixelSample, apertureSample);


			pixelColor += calculateColor(scene, sampler, aaRay);	// Calculate color for anti=aliasing ray
		}

		//Send real ray and calculate color
		pixelColor += calculateColor(scene, sampler, ray);

		//Calculate average color and return result
		return (pixelColor / (antialiaser->aaCount + 1));	// "+1" for the original ray
		
	}

	QString toString() const
	{
		return QString("Diffuse Occlusion[length=%1]").arg(m_fLength);
	}


private:
	Color3f calculateColor(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		///		!!!		BLACK = 0	WHITE = 1		!!!		///

		//Find intersection if not returns ambient color
		Intersection its;
		if(!scene->rayIntersect(ray, its))	// No intersection
		{
			return Color3f(1.f);	// White color
		}

		//Find surface normal and get costant/coefficient values
		
		//Get surface normal and light direction
		const Normal3f surfaceNormal = its.shFrame.n; //.mesh->getVertexNormals();
		Color3f pixelColor = 0.f;
		
		//Vector3f local_wo = Frame(surfaceNormal).toLocal(-ray.d);

		//Calculate color for each light source
		for(int i = 0; i < scene->getLights().size(); i++)
		{
			//Send shadow ray if it's in shadow return ambient coolor
			LightSource* currentLight = scene->getLights()[i];
			Vector3f lightDirection = currentLight->getDirection(its.p);	// Get light direction
			
			Ray3f shadowRay(its.p, lightDirection, Epsilon,  1e4f);	//Create shadow ray

			if(scene->rayIntersect(shadowRay))	// Shoots shadow ray
			{
				//It is in shadow
				pixelColor += 0.f; //+= ambientColor;
				continue;	// No need to make color calculation for this light source
			}

			//Vector3f local_wi = Frame(surfaceNormal).toLocal(lightDirection);

			//Calculate theta = N dot L
			float theta =  std::max(surfaceNormal.dot(lightDirection), 0.f); //local_wi);

			//Get Ia, Ka, Ip, and Kd

			//Calculate the color using I = IaKa + IpKd(N dot L) formulate
			BSDFQueryRecord bRec(lightDirection); //local_wi, local_wo, ESolidAngle);
			
			pixelColor += currentLight->getIntensity() * its.mesh->getBSDF()->eval(bRec) * theta;
		}

		//std::cout<<"COLOR x: "<<pixelColor.x()<<"   Y: "<<pixelColor.y()<<"   Z: "<<pixelColor.z();

		if (!pixelColor.isValid())
			return Color3f(0.0f);

		return pixelColor;
	}


	float m_fLength;
};
#pragma optimize("", on)
NORI_REGISTER_CLASS(DiffuseOcclusion, "do");
NORI_NAMESPACE_END