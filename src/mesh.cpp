/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2012 by Wenzel Jakob and Steve Marschner.

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <Eigen/Geometry>

#define NORI_TRICLIP_MAXVERTS 10

NORI_NAMESPACE_BEGIN

#pragma optimize("", off)

Mesh::Mesh() : m_vertexPositions(0), m_vertexNormals(0),
  m_vertexTexCoords(0), m_indices(0), m_vertexCount(0),
  m_triangleCount(0), m_bsdf(NULL) { }

Mesh::Mesh(const PropertyList &propList)
{
	m_bIsLight	=	propList.getBoolean("isLight", false);
	m_iId		=	propList.getInteger("id", -1);
	m_fKd		=	propList.getFloat("kd", 0.99);
}


Mesh::~Mesh() {
	//delete[] m_vertexPositions;

	if(m_vertexPositions)
	{
		delete[] m_vertexPositions;
	}

	if (m_vertexNormals)
		delete[] m_vertexNormals;
	if (m_vertexTexCoords)
		delete[] m_vertexTexCoords;
	delete[] m_indices;

	if (m_bsdf)
		delete m_bsdf;
}

void Mesh::activate() {
	/* Create a discrete distribution for sampling triangles
	   with respect to their surface area */
	m_distr.clear();
	m_distr.reserve(m_triangleCount);
	for (uint32_t i=0; i<m_triangleCount; ++i)
		m_distr.append(surfaceArea(i));
	m_distr.normalize();

	if (!m_bsdf) {
		/* If no material was assigned, instantiate a diffuse BRDF */
		m_bsdf = static_cast<BSDF *>(
			NoriObjectFactory::createInstance("diffuse", PropertyList()));
	}
}

void Mesh::samplePosition(const Point2f &_sample, Point3f &p, Normal3f &n) const {
	Point2f sample(_sample);

	/* First, sample a triangle with respect to surface area */
	size_t index = m_distr.sampleReuse(sample.x());

	/* Lookup vertex positions for the chosen triangle */
	int i0 = m_indices[3*index],
		i1 = m_indices[3*index+1],
		i2 = m_indices[3*index+2];

	const Point3f
		&p0 = m_vertexPositions[i0],
		&p1 = m_vertexPositions[i1],
		&p2 = m_vertexPositions[i2];

	/* Sample a position in barycentric coordinates */
	Point2f b = squareToUniformTriangle(sample);
	
	p = p0 * (1.0f - b.x() - b.y()) + p1 * b.x() + p2 * b.y();

	/* Also provide a normal (interpolated if vertex normals are provided) */
	if (m_vertexNormals) {
		const Normal3f 
			&n0 = m_vertexNormals[i0],
			&n1 = m_vertexNormals[i1],
			&n2 = m_vertexNormals[i2];
		n = (n0 * (1.0f - b.x() - b.y()) + n1 * b.x() + n2 * b.y()).normalized();
	} else {
		n = (p1-p0).cross(p2-p0).normalized();
	}
}

float Mesh::surfaceArea(uint32_t index) const {
	int i0 = m_indices[3*index],
		i1 = m_indices[3*index+1],
		i2 = m_indices[3*index+2];

	const Point3f
		&p0 = m_vertexPositions[i0],
		&p1 = m_vertexPositions[i1],
		&p2 = m_vertexPositions[i2];
	
	return 0.5f * Vector3f((p1-p0).cross(p2-p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
	int i0 = m_indices[3*index],
		i1 = m_indices[3*index+1],
		i2 = m_indices[3*index+2];

	const Point3f
		&p0 = m_vertexPositions[i0],
		&p1 = m_vertexPositions[i1],
		&p2 = m_vertexPositions[i2];

	/* find vectors for two edges sharing v[0] */
	Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

	/* begin calculating determinant - also used to calculate U parameter */
	Vector3f pvec = ray.d.cross(edge2);

	/* if determinant is near zero, ray lies in plane of triangle */
	float det = edge1.dot(pvec);

	if (det > -1e-8f && det < 1e-8f)
		return false;
	float inv_det = 1.0f / det;

	/* calculate distance from v[0] to ray origin */
	Vector3f tvec = ray.o - p0;

	/* calculate U parameter and test bounds */
	u = tvec.dot(pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	/* prepare to test V parameter */
	Vector3f qvec = tvec.cross(edge1);

	/* calculate V parameter and test bounds */
	v = ray.d.dot(qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	/* ray intersects triangle -> compute t */
	t = edge2.dot(qvec) * inv_det;

	return true;
}
	
BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
	BoundingBox3f result(m_vertexPositions[m_indices[3*index]]);
	result.expandBy(m_vertexPositions[m_indices[3*index+1]]);
	result.expandBy(m_vertexPositions[m_indices[3*index+2]]);
	return result;
}

/// Internally used by getClippedBoundingBox()
static int sutherlandHodgman(Point3d *input, int inCount, Point3d *output, int axis, 
		double splitPos, bool isMinimum) {
	if (inCount < 3)
		return 0;

	Point3d cur       = input[0];
	double sign       = isMinimum ? 1.0f : -1.0f;
	double distance   = sign * (cur[axis] - splitPos);
	bool  curIsInside = (distance >= 0);
	int   outCount    = 0;

	for (int i=0; i<inCount; ++i) {
		int nextIdx = i+1;
		if (nextIdx == inCount)
			nextIdx = 0;
		const Point3d &next = input[nextIdx];
		distance = sign * (next[axis] - splitPos);
		bool nextIsInside = (distance >= 0);

		if (curIsInside && nextIsInside) {
			/* Both this and the next vertex are inside, add to the list */
			assert(outCount + 1 < NORI_TRICLIP_MAXVERTS);
			output[outCount++] = next;
		} else if (curIsInside && !nextIsInside) {
			/* Going outside -- add the intersection */
			double t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			assert(outCount + 1 < NORI_TRICLIP_MAXVERTS);
			Point3d p = cur + (next - cur) * t;
			p[axis] = splitPos; // Avoid roundoff errors
			output[outCount++] = p;
		} else if (!curIsInside && nextIsInside) {
			/* Coming back inside -- add the intersection + next vertex */
			double t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			assert(outCount + 2 < NORI_TRICLIP_MAXVERTS);
			Point3d p = cur + (next - cur) * t;
			p[axis] = splitPos; // Avoid roundoff errors
			output[outCount++] = p;
			output[outCount++] = next;
		} else {
			/* Entirely outside - do not add anything */
		}
		cur = next;
		curIsInside = nextIsInside;
	}
	return outCount;
}

BoundingBox3f Mesh::getClippedBoundingBox(uint32_t index, const BoundingBox3f &bbox) const {
	/* Reserve room for some additional vertices */
	Point3d vertices1[NORI_TRICLIP_MAXVERTS], vertices2[NORI_TRICLIP_MAXVERTS];
	int nVertices = 3;

	/* The kd-tree code will frequently call this function with
	   almost-collapsed bounding boxes. It's extremely important not to introduce
	   errors in such cases, otherwise the resulting tree will incorrectly
	   remove triangles from the associated nodes. Hence, do the
	   following computation in double precision! */
	for (int i=0; i<3; ++i) 
		vertices1[i] = m_vertexPositions[m_indices[3*index+i]].cast<double>();

	for (int axis=0; axis<3; ++axis) {
		nVertices = sutherlandHodgman(vertices1, nVertices, vertices2, axis, bbox.min[axis], true);
		nVertices = sutherlandHodgman(vertices2, nVertices, vertices1, axis, bbox.max[axis], false);
	}

	BoundingBox3f result;
	for (int i=0; i<nVertices; ++i) 
		result.expandBy(vertices1[i].cast<float>());
	result.clip(bbox);
	return result;
}

void Mesh::addChild(NoriObject *obj) {
	switch (obj->getClassType()) {
		case EBSDF:
			if (m_bsdf)
				throw NoriException("Mesh: tried to register multiple BSDF instances!");
			m_bsdf = static_cast<BSDF *>(obj);
			break;

		default:
			throw NoriException(QString("Mesh::addChild(<%1>) is not supported!").arg(
				classTypeName(obj->getClassType())));
	}
}

QString Mesh::toString() const {
	return QString(
		"Mesh[\n"
		"  name = \"%1\",\n"
		"  vertexCount = %2,\n"
		"  triangleCount = %3,\n"
		"  bsdf = %4\n"
		"]")
	.arg(m_name)
	.arg(m_vertexCount)
	.arg(m_triangleCount)
	.arg(indent(m_bsdf->toString()));
}

QString Intersection::toString() const {
	if (!mesh)
		return "Intersection[invalid]";

	return QString(
		"Intersection[\n"
		"  p = %1,\n"
		"  t = %2,\n"
		"  uv = %3,\n"
		"  shFrame = %4,\n"
		"  geoFrame = %5,\n"
		"  mesh = %6\n"
		"]")
	.arg(p.toString())
	.arg(t)
	.arg(uv.toString())
	.arg(indent(shFrame.toString()))
	.arg(indent(geoFrame.toString()))
	.arg(indent(mesh->toString()));
}

#pragma optimize("", on)

NORI_NAMESPACE_END
