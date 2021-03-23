#include "pzgmesh.h"//geometrical mesh
#include "TPZVTKGeoMesh.h"//for printing the mesh in .vtk format
//special maps
#include "TPZQuadTorus.h"
#include "TPZTriangleTorus.h"
#include "tpzellipse3d.h"
#include "tpzquadraticline.h"
#include "tpzquadraticpyramid.h"
#include "tpzquadratictetra.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "TPZWavyLine.h"
#include "tpzarc3d.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"


/*
This function will first create a mesh with just one element
of the type TGeo. Then, this mesh is refined nref times
for better visualization.
*/
template <class TGeo>
void CreateSampleElement(const int nref);

int main(int argc, char *argv[]) {
    constexpr int nRef = 4;

    /*The first two classes create an element based on a torus
      whose axis is aligned with the z-direction.*/
    
    //Creating a (incomplete) torus from a quadrilateral
    CreateSampleElement<pzgeom::TPZQuadTorus>(nRef);
    //Creating a (incomplete) torus from a triangle.
    CreateSampleElement<pzgeom::TPZTriangleTorus>(nRef);

    //Creating an element that maps a lign segment to an ellipse arc
    CreateSampleElement<pzgeom::TPZEllipse3D>(nRef);

    //Creating one half of a spherical shell from a quadrilateral
    CreateSampleElement<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad>>(nRef);
    //Creating one quarter of a spherical shell from a triangle
    CreateSampleElement<pzgeom::TPZTriangleSphere<pzgeom::TPZGeoTriangle>>(nRef);
    //Mapping a line segment to a cosine
    CreateSampleElement<pzgeom::TPZWavyLine>(nRef);
    //Mapping a line segment to an arc.
    
    /*This element is really useful interfaces in 2D problems
      with curved geometries.*/
    CreateSampleElement<pzgeom::TPZArc3D>(nRef);

    /*The following are elements obtained using through
      quadratic lagrangian functions*/
    CreateSampleElement<pzgeom::TPZQuadraticLine>(nRef);
    CreateSampleElement<pzgeom::TPZQuadraticTrig>(nRef);
    CreateSampleElement<pzgeom::TPZQuadraticQuad>(nRef);
    CreateSampleElement<pzgeom::TPZQuadraticTetra>(nRef);
    CreateSampleElement<pzgeom::TPZQuadraticCube>(nRef);
    CreateSampleElement<pzgeom::TPZQuadraticPyramid>(nRef);
    CreateSampleElement<pzgeom::TPZQuadraticPrism>(nRef);
}




template <class TGeo>
void CreateSampleElement(const int nref){
    std::string elName = TGeo::TypeName();

    std::cout<<"Starting mesh refinement of element type: "<<elName<<std::endl;

    //material identifier
    constexpr int matId = 1;
    /*
      Geometric coordinates in NeoPZ are always assumed to be in a
      three dimensional space.

      The function of the parameters lowerCorner and size will differ
      from element to element.

      TPZManVector<TYPE,N> are vectors with bound checking and static
      allocation (with dinamic allocation being possible if needed)*/
    TPZManVector<REAL,3> lowerCorner(3,0);
    TPZManVector<REAL,3> size(3,1);
    //Geometrical mesh
    TPZGeoMesh gmesh;
    TGeo::InsertExampleElement(gmesh, matId, lowerCorner, size);
    /*
      Computes the neighbouring relationship in the mesh.
      Necessary for the refinement process.
     */
    gmesh.BuildConnectivity();

    /*In the following loop, each element that has not been
     divided (this check is performed through
     TPZGeoEl::NSubElements()) will be refined using its default
     refinement pattern.*/
    {
        TPZVec<TPZGeoEl *> sons;
        for(int iDiv = 0; iDiv < nref; iDiv++){
            const int nel = gmesh.NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh.Element(iel);
                if (geo->NSubElements() == 0) {
                    geo->Divide(sons);
                }
            }
        }
    }

    //outputs the mesh in .vtk format
    std::ofstream fileName(elName+".vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&gmesh, fileName);
}