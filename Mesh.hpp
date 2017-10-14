#ifndef __MESH_HPP__

#include <QtGui>
#include <QtOpenGL>
#include <algorithm>
#include <complex>

#include <iostream>
using namespace std;

struct Mesh_Face {
    Mesh_Face() {
        vert[0] = vert[1] = vert[2] = -1;
    }

    Mesh_Face(long v0, long v1, long v2) {
        vert[0] = v0; vert[1] = v1; vert[2] = v2;
    }
    long vert[3]; // indices (in the vertex array) of all vertices (mesh_vertex)

    long adj[3]; // neigbor faces
    QVector3D FN;
};

struct Mesh_Vertex {
    QVector3D me;
    QVector3D normal;
    double avg_edge;
    vector<long> neighbor;
    vector<long> face_neighbor;
    QVector3D gaussian;
    QVector3D laplacian;
    int count;
};

struct Mesh {
    vector<QVector3D> vertices; // List of shared verticies.
    vector<Mesh_Face> faces; // Mesh faces.
    QOpenGLBuffer vertexBuffer, baryBuffer;
    vector<Mesh_Vertex> vertices_list;

    void avg_length();
    void process_all();
    bool load_obj(QString filename);
    void storeVBO();
    void recenter();
    void add_face(const vector<int> &cur_vert);
    void process_example();
    void inflate(double factor);
    void compute_allNormals();
    void randomNoise(double factor);
    void gaussian();
    void smooth();
    void sharpen();
    void laplacian();
    void set_up_vertices_list();
    void centerVerticesTangentially();
    long intersect(long v1, long v2, long curr);
    void find_adjFace();
    void splitFaces();
    void splitLongEdges();
    long findSingleFace(long e, long uncommon, long f);
    void splitWork(double length, double threshold, long f, long f2,
                               long e1, long e2);
    long findUncommonPt(long f, long e1, long e2);
    long insert_midPT(long v1, long v2);
    long insert_midPT(long v1, long v2, vector<QVector3D> &buffer, vector<QVector3D> &ori);
    void calculate_evenBeta(long index, vector<QVector3D> &buffer,
                                       vector<QVector3D> &curr);
    void loopSubdivision();

};

#endif // __MESH_HPP__


