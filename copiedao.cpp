#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class CopiedAO : public Integrator
{
public:
	CopiedAO(const PropertyList &propList) {
		/* Ray length of the ambient occlusion queries;
		   expressed relative to the scene size */
		m_length = propList.getFloat("length", 0.1f);
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, int indeX = 0, int indexY = 0) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		/* Sample a cosine-weighted direction from the hemisphere (local coordinates) */
		Vector3f d = squareToCosineHemisphere(sampler->next2D());

		/* Use the shading frame at "its" to convert it to world coordinates */
		d = its.toWorld(d);

		/* Determine the length of the "shadow ray" based on the scene size
		   and the configuration options */
		float length = m_length * scene->getBoundingBox().getExtents().norm();

		/* Create a new outgoing ray having extents (epsilon, length) */
		Ray3f shadowRay(its.p, d, Epsilon, length);

		/* Perform an occlusion test and return one or zero depending on the result */
		return Color3f(scene->rayIntersect(shadowRay) ? 0.0f : 1.0f);
	}

	QString toString() const {
		return QString("CopiedAmbientOcclusion[length=%1]").arg(m_length);
	}
private:
	float m_length;
};

NORI_REGISTER_CLASS(CopiedAO, "copiedao");
NORI_NAMESPACE_END