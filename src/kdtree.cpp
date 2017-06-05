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

#include <nori/kdtree.h>
#include <Eigen/Geometry>

#include "..\sphere.h"

NORI_NAMESPACE_BEGIN
#pragma optimize( "", on )
KDTree::KDTree() : m_primitiveCount(0) {
	m_sizeMap.push_back(0);
}

KDTree::~KDTree() {
	for (size_t i=0; i<m_meshes.size(); ++i)
		delete m_meshes[i];
}

void KDTree::build() {
	SizeType primCount = getPrimitiveCount();
	cout << "Constructing a SAH kd-tree (" << primCount << " triangles, "
		 << getCoreCount() << " threads) .." << endl;
	Parent::buildInternal();
}

void KDTree::addMesh(Mesh *mesh) {
	m_primitiveCount += mesh->getTriangleCount();
	m_meshes.push_back(mesh);
	m_sizeMap.push_back(m_sizeMap.back() + mesh->getTriangleCount());
}

bool KDTree::rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const {
	/// KD-tree traversal stack
	struct {
		/* Pointer to the far child */
		const KDNode * __restrict node;
		/* Distance traveled along the ray (entry or exit) */
		float t;
		/* Previous stack item */
		uint32_t prev;
		/* Associated point */
		Point3f p;
	} stack[NORI_KD_MAXDEPTH];

	its.t = std::numeric_limits<float>::infinity();

	/* Use an adaptive ray epsilon */
	float mint = ray.mint, maxt = ray.maxt;
	if (mint == Epsilon) 
		mint = std::max(mint, mint * ray.o.array().abs().maxCoeff());

	float bboxMinT, bboxMaxT;
	if (!m_bbox.rayIntersect(ray, bboxMinT, bboxMaxT))
		return false;

	mint = std::max(mint, bboxMinT);
	maxt = std::min(maxt, bboxMaxT);

	if (maxt < mint)
		return false;

	/* Set up the entry point */
	uint32_t enPt = 0;
	stack[enPt].t = mint;
	stack[enPt].p = ray(mint);

	/* Set up the exit point */
	uint32_t exPt = 1;
	stack[exPt].t = maxt;
	stack[exPt].p = ray(maxt);
	stack[exPt].node = NULL;

	bool foundIntersection = false;
	uint32_t foundPrimIndex = 0;
	const KDNode * __restrict currNode = m_nodes;
	while (currNode != NULL) {
		while (EXPECT_TAKEN(!currNode->isLeaf())) {
			const float splitVal = (float) currNode->getSplit();
			const int axis = currNode->getAxis();
			const KDNode * __restrict farChild;

			if (stack[enPt].p[axis] <= splitVal) {
				if (stack[exPt].p[axis] <= splitVal) {
					/* Cases N1, N2, N3, P5, Z2 and Z3 (see thesis) */
					currNode = currNode->getLeft();
					continue;
				}

				/* Typo in Havran's thesis:
				   (it specifies "stack[exPt].p == splitVal", which
				    is clearly incorrect) */
				if (stack[enPt].p[axis] == splitVal) {
					/* Case Z1 */
					currNode = currNode->getRight();
					continue;
				}

				/* Case N4 */
				currNode = currNode->getLeft();
				farChild = currNode + 1; // getRight()
			} else { /* stack[enPt].p[axis] > splitVal */
				if (splitVal < stack[exPt].p[axis]) {
					/* Cases P1, P2, P3 and N5 */
					currNode = currNode->getRight();
					continue;
				}
				/* Case P4 */
				farChild = currNode->getLeft();
				currNode = farChild + 1; // getRight()
			}

			/* Cases P4 and N4 -- calculate the distance to the split plane */
			float distToSplit = (splitVal - ray.o[axis]) * ray.dRcp[axis];

			/* Set up a new exit point */
			const uint32_t tmp = exPt++;
			if (exPt == enPt) /* Do not overwrite the entry point */
				++exPt;

			stack[exPt].prev = tmp;
			stack[exPt].t = distToSplit;
			stack[exPt].node = farChild;
			stack[exPt].p = ray(distToSplit);
			stack[exPt].p[axis] = splitVal;
		}

		/* Reached a leaf node */
		for (IndexType entry=currNode->getPrimStart(),
				last = currNode->getPrimEnd(); entry != last; entry++) {
			IndexType primIndex = m_indices[entry];
			IndexType meshIndex = findMesh(primIndex);
			const Mesh *mesh = m_meshes[meshIndex];

			float u, v, t;
			bool success = mesh->rayIntersect(primIndex, ray, u, v, t);

		//	std::cout<<"KD: SUCCESS -> "<<success<<std::endl;

			//Normal3f n;
			//mesh->myRayIntersect(prim..., n);

			if (success && t >= mint && t <= maxt) {
				if (shadowRay)
					return true;
				maxt = t;
				its.t = t;
				its.uv = Point2f(u, v);
				its.mesh = mesh;
				foundPrimIndex = primIndex;
				foundIntersection = true;
			}
		}

		if (stack[exPt].t > maxt) 
			break;

		/* Pop from the stack and advance to the next node on the interval */
		enPt = exPt;
		currNode = stack[exPt].node;
		exPt = stack[enPt].prev;
	}

	if (foundIntersection && !shadowRay) {

		const Mesh *mesh = its.mesh;

		if( GeometricType == mesh->getMeshType() )
		{
			//Calculate and set intersection point
			/*
			For the calculation, Ri = [xi, yi, zi] = [x0 + xd * ti ,  y0 + yd * ti,  z0 + zd * ti]
			formula of Siggraph tutorial is used. Where, Ri is intersection point,
			where X0, Y0, Z0 are origin of the ray, Xd, Yd, and Zd are direction of the Ray and ti is t.
			*/
						//X			   //Y			  //Z
			float Xd = ray.d.x();		float Yd = ray.d.y();		float Zd = ray.d.z();
			float Xo = ray.o.x();		float Yo = ray.o.y();		float Zo = ray.o.z();

			//its.p is the intersection point/position
			its.p.x() = Xo + Xd * its.t;
			its.p.y() = Yo + Yd * its.t;
			its.p.z() = Zo + Zd * its.t;

			// SURFACE NORMAL CALCULATION

			//Calculate intersection position
			
			//Point3f Ri((Xo + Xd * its.t), (Yo + Yd * its.t), (Zo + Zd * its.t)); 

			//Calculate Surface Normal
			/*
			For the calculation, SN = [(xi - xc)/Sr,   (yi - yc)/Sr,   (zi - zc)/Sr]
			formula from Siggraph tutorial is used. Where, SN is the surface normal,
			Xi, Yi, Zi are intersection point coordinates, Xc, Yc, Zc are coordinates
			of the sphere's center.
			*/

			const Sphere *currentSphere = (Sphere*)its.mesh;
			const float radius = currentSphere->getRadius();
			const Point3f center = currentSphere->getCenter();

			float invRadius = 1.f / radius;	// Calculate inverse radius. Because * is faster than /

			//Calculate normal
			float Nx = (its.p.x() - center.x()) * invRadius;
			float Ny = (its.p.y() - center.y()) * invRadius;
			float Nz = (its.p.z() - center.z()) * invRadius;

			//Assing normal to Intersection object
			its.shFrame.n = Normal3f(Nx, Ny, Nz);
			its.shFrame.n.normalize();

			//Calculate UV Mapping values

			/*
			UV Mapping calculation formulae are written with referenced Wikipedia Page
			Link:	http://en.wikipedia.org/wiki/UV_mapping
			*/

			//Calculate Vector d, which is unit vector from intersection point to the sphere's origin

			// Intersection point Ri can be calculated using the equation below (Siggraph article):
			// Ri = [xi, yi, zi] = [x0 + xd * ti ,  y0 + yd * ti,  z0 + zd * ti]

			//Calculate Vector d
			Vector3f d = its.p - center;
			d.normalize();	// Normalize Vector d

			//Calculate U value
			float u = 0.5 + (std::atan2(d.z(), d.x()) / (2 * M_PI));

			//Calculte V value
			float v = 0.5 - (std::asin(d.y()) / M_PI);

			//Assing UV values to Intersection object
			its.uv = Point2f(u, v);

			return foundIntersection;
		}


		/* Find the barycentric coordinates */
		Vector3f bary;
		bary << 1-its.uv.sum(), its.uv;

		/* Look up the vertex indices */
		//const Mesh *mesh = its.mesh;
		const uint32_t *indices = mesh->getIndices(),
				  idx0 = indices[3*foundPrimIndex+0],
				  idx1 = indices[3*foundPrimIndex+1],
				  idx2 = indices[3*foundPrimIndex+2];

		const Point3f  *positions = mesh->getVertexPositions();
		const Normal3f *normals   = mesh->getVertexNormals();
		const Point2f  *texCoords = mesh->getVertexTexCoords();

		Point3f p0 = positions[idx0],
			p1 = positions[idx1],
			p2 = positions[idx2];

		/* Compute the intersection positon accurately 
		   using barycentric coordinates */
		its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

		/* Compute proper texture coordinates if provided by the mesh */
		if (texCoords) 
			its.uv = bary.x() * texCoords[idx0] +
				bary.y() * texCoords[idx1] +
				bary.z() * texCoords[idx2];

		/* Compute the geometry frame */
		its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

		if (normals) {
			/* Compute the shading frame. Note that for simplicity,
			   the current implementation doesn't attempt to provide
			   tangents that are continuous across the surface. That
			   means that this code will need to be modified to be able
			   use anisotropic BRDFs, which need tangent continuity */

			its.shFrame = Frame(
				(bary.x() * normals[idx0] +
				 bary.y() * normals[idx1] +
				 bary.z() * normals[idx2]).normalized());
		} else {
			its.shFrame = its.geoFrame;
		}
	}

	//std::cout<<"KD: foundIntersection -> "<<foundIntersection<<std::endl;

	return foundIntersection;
}
#pragma optimize( "", off )
NORI_NAMESPACE_END
