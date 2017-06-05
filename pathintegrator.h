/*
	DEVELOPER:	MCO

	Path tracing integrator
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <nori\bsdf.h>	// For using BSDFQueryRecord object in the case of intersecting with light source

NORI_NAMESPACE_BEGIN

class PathIntegrator : public Integrator
{
public:

	virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, int indexX = 0, int indexY = 0) const = 0;

	QString toString() const{	return QString("Path Tracer");	}

protected:

	virtual Color3f tracePath(const Scene *scene, Sampler *sampler, const Intersection firstIntersection, const Vector3f firstRayDirection) const = 0;
};

NORI_NAMESPACE_END