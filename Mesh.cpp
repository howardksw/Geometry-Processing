#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <complex>

#include "Mesh.hpp"

using namespace std;

void Mesh::process_all() {
    avg_length();
    compute_allNormals();
    gaussian();
    laplacian();
    find_adjFace();
}

void Mesh::set_up_vertices_list(){
     for(uint i = 0; i < vertices.size(); i++) {
         Mesh_Vertex temp;
         temp.me = QVector3D(vertices[i]);
         vertices_list.push_back(temp);
     }
}

void Mesh::avg_length(){
    // Setting up the new data structure for Mesh_Vertex
    for(uint i = 0; i < vertices.size(); i++) {
        Mesh_Vertex temp;
        temp.me = QVector3D(vertices[i]);

        // setting up the Mesh_Vertex data structure
        for(long f = 0; f < (long)faces.size(); f++) {
            // adding vertex and face neighbors
            if (faces[f].vert[0] == i) {
                temp.neighbor.push_back(faces[f].vert[1]);
                temp.neighbor.push_back(faces[f].vert[2]);
                temp.face_neighbor.push_back(f);
            }
            else if (faces[f].vert[1] == i) {
                temp.neighbor.push_back(faces[f].vert[0]);
                temp.neighbor.push_back(faces[f].vert[2]);
                temp.face_neighbor.push_back(f);
            }
            else if (faces[f].vert[2] == i) {
                temp.neighbor.push_back(faces[f].vert[1]);
                temp.neighbor.push_back(faces[f].vert[0]);
                temp.face_neighbor.push_back(f);
            }
            // adding neighbor faces into a data structure from a vertex
        }

        // calculate edge
        if(temp.neighbor.size() != 0) {
            sort(temp.neighbor.begin(),temp.neighbor.end());

            vector<long>::iterator it;
            it = unique(temp.neighbor.begin(), temp.neighbor.end());
            temp.neighbor.resize(distance(temp.neighbor.begin(),it));

            QVector3D v0 = vertices[i];
            double sum_edge = 0.0;
            long index = 0;
            for (long ele:temp.neighbor) {
                QVector3D v1 = vertices[ele] - v0;
                index++;
                sum_edge += qSqrt(qPow(v1.x(),2) + qPow(v1.y(),2) + qPow(v1.z(),2));
            }
            temp.avg_edge = sum_edge / (long)index;
            temp.normal = QVector3D(0,0,0);
        }
        // update
        vertices_list[i] = temp;
    }
}

void Mesh::compute_allNormals() {
    // Find all the faces' normals
    for(long i = 0; i < (long) faces.size(); i++) {
        // Finding 3 vertices on a triangle
        QVector3D v0 = vertices[faces[i].vert[0]];
        QVector3D v1 = vertices[faces[i].vert[1]];
        QVector3D v2 = vertices[faces[i].vert[2]];

        // find the 2 vectors and do cross product
        QVector3D product = QVector3D::crossProduct((v1 - v0), (v2 - v0));

        faces[i].FN = product;

        // Compute all the vertices' normals

       vertices_list[faces[i].vert[0]].normal += faces[i].FN;
       vertices_list[faces[i].vert[1]].normal += faces[i].FN;
       vertices_list[faces[i].vert[2]].normal += faces[i].FN;
    }

    //normalize the all vertices
    for(uint i = 0; i < vertices_list.size(); i++) {
        vertices_list[i].normal = vertices_list[i].normal.normalized();
    }
}


void Mesh::inflate(double factor) {
    // iterate through the vertices vector
    for(uint i = 0; i < vertices.size(); i++) {
        // get the new distance by multiplying the factor and the average edge
        // on the current vertex
        double distance = vertices_list[i].avg_edge * factor;
        // get the new vector by multiplying the normal vector
        QVector3D new_vector = vertices_list[i].normal;

        // add the new_vector with the original vertex
        vertices[i] += (distance * new_vector);
    }
}

void Mesh::randomNoise(double factor) {
    for (long i = 0; i < (long)  vertices_list.size(); i++) {
        // find the maximum range
        double maximum = factor * vertices_list[i].avg_edge;
        // calculate random number
        double random = rand() % int((maximum + 1) * 1000);
        random = random / 1000.0;

        // find the random increasing rate on xyz
        double rand_x = rand() % (1001);
        double rand_y = rand() % (1001);
        double rand_z = rand() % (1001);

        // update with the new positions
        vertices[i] += (QVector3D(rand_x,rand_y,rand_z) / 10000.0) * random;
    }
}

void Mesh::gaussian() {
    double pi = 3.14159265;
    // Iterate through all the vertices
    for (long i = 0; i < (long)  vertices_list.size(); i++) {
        double sigma = vertices_list[i].avg_edge;
        double sum_weight = 0.0;
        QVector3D sum = QVector3D(0,0,0);

        // Iterate through all the neighbors attached to this specific vertex
        for(long n:vertices_list[i].neighbor) {
            double x = vertices[i].x() - vertices[n].x();
            double y = vertices[i].y() - vertices[n].y();
            double z = vertices[i].z() - vertices[n].z();

            // gaussian function
            double e = exp(-(qPow(x,2)+qPow(y,2) + qPow(z,2))/(2*qPow(sigma,2)));
            double constant = 1/(qSqrt((qPow(2*pi,3)))*qPow(sigma,3));
            double weight = e * constant;

            // add up the weighted results
            sum += vertices[n] * weight;
            sum_weight += weight;
        }

        //for itself
        double constant = 1 / (qSqrt((qPow(2*pi,3)))*qPow(sigma,3));
        double weight = 1 * constant;
        sum_weight += weight;
        sum += vertices[i] * weight;

        // update
        QVector3D result = sum / sum_weight;
        vertices_list[i].gaussian = result;
    }
}

void Mesh::laplacian() {
    // Iterate through all the vertices
    for (long i = 0; i < (long)  vertices_list.size(); i++) {

        long count = 0;
        QVector3D sum = QVector3D(0,0,0);
        // Iterate through all the neighbors attached to this specific vertex
        for(long n:vertices_list[i].neighbor) {
            count++;
            sum += vertices[n];
        }
        //for itself
        count++;
        sum += vertices[i];

        QVector3D result = sum / count;
        vertices_list[i].laplacian = result;
    }
}

// intersect function to find the index of the intersected face of two vertices
long Mesh::intersect(long v1, long v2, long curr) {
    for (long i = 0; i < (long) vertices_list[v1].face_neighbor.size(); i++) {
        for(long j = 0; j < (long) vertices_list[v2].face_neighbor.size(); j++) {
            // when both neighbor faces are the same and not equal to the current face
            if (vertices_list[v1].face_neighbor[i] == vertices_list[v2].face_neighbor[j]
                    && curr != vertices_list[v1].face_neighbor[i]) {
                return (long) vertices_list[v1].face_neighbor[i];
            }
        }
    }
    return -1;
}

// build adjacency data structure
void Mesh::find_adjFace() {
    for (long i = 0; i < (long) faces.size(); i++) {
        faces[i].adj[0] = intersect(faces[i].vert[0], faces[i].vert[1], i);
        faces[i].adj[1] = intersect(faces[i].vert[1], faces[i].vert[2], i);
        faces[i].adj[2] = intersect(faces[i].vert[0], faces[i].vert[2], i);
    }
}

// smooth the mesh with gaussian
void Mesh::smooth() {
    for (long i = 0; i < (long) vertices_list.size(); i++) {
        vertices[i] = vertices_list[i].gaussian;
    }
}

void Mesh::sharpen() {}

// Using laplacian function to update to a new mesh
void Mesh::centerVerticesTangentially() {
    for (long i = 0 ; i < (long) vertices_list.size(); i++) {
        vertices[i] = vertices_list[i].laplacian;
    }
}

long Mesh::insert_midPT(long v1, long v2) {
    // find the mid point of an edge
    QVector3D mid = (vertices[v1] + vertices[v2]) / 2;

    // check if the mid point is contained in the vertices list
    for(uint i = 0; i < vertices.size(); i++) {
        // if yes, return the current index
        if (mid == vertices[i]) {
            return i;
        }
    }
    vertices.push_back(mid);
    return (long) vertices.size() - 1;
}

long Mesh::insert_midPT(long v1, long v2, vector<QVector3D> &buffer, vector<QVector3D> &ori) {
    // find the mid point of an edge
    QVector3D mid = (buffer[v1] + buffer[v2]) / 2;

    // check if the mid point is contained in the vertices list
    for(uint i = 0; i < ori.size(); i++) {
        // if yes, return the current index
        if (mid == ori[i]) {
            return i;
        }
    }
    ori.push_back(mid);
    return (long) vertices.size() - 1;
}

void Mesh::splitFaces() {
    long count = faces.size();
    for (long i = 0; i < count; i++) {
        long v0 = faces[i].vert[0];
        long v1 = faces[i].vert[1];
        long v2 = faces[i].vert[2];
        // find three new mid points on the edges
        long mid1 = insert_midPT(v0,v1);
        long mid2 = insert_midPT(v0,v2);
        long mid3 = insert_midPT(v1,v2);

        // first triangle by creating a new one
        Mesh_Face f1 = Mesh_Face(v0, mid1, mid2);

        // second triangle by creating a new one
        Mesh_Face f2 = Mesh_Face(mid1, v1, mid3);

        // third triangle by creating a new one
        Mesh_Face f3 = Mesh_Face(mid2, mid3, v2);

        // fourth triangle by changing the original triangle
        faces[i].vert[0] = mid1;  // mid1
        faces[i].vert[1] = mid2;
        faces[i].vert[2] = mid3;

        // add the new faces
        faces.push_back(f1);
        faces.push_back(f2);
        faces.push_back(f3);
    }
}

long Mesh::findUncommonPt(long f, long e1, long e2) {
    // find the uncommon point from two neighbor triangles
    for(long i = 0; i < 3; i++) {
        if ((faces[f].vert[i] != e1 && faces[f].vert[i] != e2)) {
            return faces[f].vert[i];
        }
    }
    return -1;
}

// find an adjacency face based on a edge and a face
long Mesh::findSingleFace(long e, long uncommon, long f) {
    long adj;
    for (long j = 0; j < 3; j++) {
        if ((faces[f].vert[0] == e && faces[f].vert[1] == uncommon) ||
             (faces[f].vert[1] == e && faces[f].vert[0] == uncommon)) {
            adj = faces[f].adj[0];

        }else if ((faces[f].vert[1] == e && faces[f].vert[2] == uncommon) ||
                  (faces[f].vert[2] == e && faces[f].vert[1] == uncommon)) {
            adj = faces[f].adj[1];

        }else if ((faces[f].vert[0] == e && faces[f].vert[2] == uncommon) ||
                  (faces[f].vert[2] == e && faces[f].vert[0] == uncommon)) {
            adj = faces[f].adj[2];
        }
    }
    return adj;
}

// splitWork is called by splitLongEdges
// doing the actual splitting operation on two faces
void Mesh::splitWork(double length, double threshold, long f, long f2, long e1, long e2) {
    if (length > threshold) {
        // find the mid point
        long pos = insert_midPT(e1, e2);
        // find the uncommon points on two neighbor faces
        long f1_uncommon = findUncommonPt(f, e1, e2);
        long f2_uncommon = findUncommonPt(f2, e1, e2);

        // making new faces
        Mesh_Face new_face = Mesh_Face(e1,pos,f1_uncommon);
        Mesh_Face new_face2 = Mesh_Face(e1,pos,f2_uncommon);

        // get all the original adjacency faces
        long adj1 = findSingleFace(e1,f1_uncommon,f);
        long adj2 = findSingleFace(e2,f1_uncommon,f);
        long adj3 = findSingleFace(e1,f2_uncommon,f2);
        long adj4 = findSingleFace(e2,f2_uncommon,f2);

        // updating the old one as the splited one
        faces[f].vert[0] = e2;
        faces[f].vert[1] = pos;
        faces[f].vert[2] = f1_uncommon;

        faces[f2].vert[0] = e2;
        faces[f2].vert[1] = pos;
        faces[f2].vert[2] = f2_uncommon;

        // adding faces to the mesh
        faces.push_back(new_face);
        faces.push_back(new_face2);

        // position of the two new faces
        long face_idx1 = faces.size() - 2;
        long face_idx2 = faces.size() - 1;

        // update adjacency faces
        new_face.adj[0] = face_idx2;
        new_face.adj[1] = f;
        new_face.adj[2] = adj1;

        new_face2.adj[0] = face_idx1;
        new_face2.adj[1] = f2;
        new_face2.adj[2] = adj3;

        faces[f].adj[0] = f2;
        faces[f].adj[1] = face_idx1;
        faces[f].adj[2] = adj2;

        faces[f2].adj[0] = f;
        faces[f2].adj[1] = face_idx2;
        faces[f2].adj[2] = adj4;
    }
}

void Mesh::splitLongEdges() {
    // firstly initialize the count_full vector for checking the neighbor vertices
    vector<int> count_full;
    if (count_full.size() == 0) {
        for(long i = 0; i < (long) vertices.size(); i++) {
            count_full.push_back(0);
        }
    }else {
        for(long i = 0; i < (long) vertices.size(); i++) {
            count_full[i] = 0;
        }
    }

    // now starting to calculate the average edge length
    double sum_edge = 0.0;
    int num_edge = 0;
    for(long i = 0; i < (long) vertices.size(); i++) {

        // iterate through the neighbor lists
        for(long j = 0; j < (long) vertices_list[i].neighbor.size(); j++) {

            long pos = vertices_list[i].neighbor[j];
            // check if the current vertices has been calculated all edges with its neighbors
            if (count_full[pos] < (int) vertices_list[pos].neighbor.size()) {
                // calculate the edge and add it to sum
                QVector3D v = vertices[i] - vertices[pos];
                sum_edge += qSqrt(qPow(v.x(),2) + qPow(v.y(),2) + qPow(v.z(),2));
                num_edge++;
            }
            count_full[i] = count_full[i] + 1;
        }
    }
    // average edges value
    double avg_edge = (sum_edge/num_edge) * 4/3;

//    long size = faces.size();
    for (long i = 0; i < (long) faces.size(); i++) {
        // Iterate through all the faces and check all the edges
        // see if they are longer than the 4/3 of the average edges value
        long v0 = faces[i].vert[0];
        long v1 = faces[i].vert[1];
        long v2 = faces[i].vert[2];

        // first edge
        QVector3D v = vertices[v0] - vertices[v1];
        double length = qSqrt(qPow(v.x(),2) + qPow(v.y(),2) + qPow(v.z(),2));

        long adj_face = faces[i].adj[0];
        // if having an adjacency face, then split
        if (adj_face != -1) {
            splitWork(length, avg_edge, i, adj_face, v0, v1);
        }

        // second edge
        v = vertices[v1] - vertices[v2];
        length = qSqrt(qPow(v.x(),2) + qPow(v.y(),2) + qPow(v.z(),2));

        adj_face = faces[i].adj[1];
        // if having an adjacency face, then split
        if (adj_face != -1) {
            splitWork(length, avg_edge, i, adj_face, v1, v2);
        }

        // third edge
        v = vertices[v0] - vertices[v2];
        length = qSqrt(qPow(v.x(),2) + qPow(v.y(),2) + qPow(v.z(),2));

        adj_face = faces[i].adj[2];
        // if having an adjacency face, then split
        if (adj_face != -1) {
            splitWork(length, avg_edge, i, adj_face, v0, v2);
        }
    }
}

// compute even vertices with beta value
void Mesh::calculate_evenBeta(long index, vector<QVector3D> &buffer,
                                   vector<QVector3D> &curr) {
    long n = vertices_list[index].neighbor.size();
    double beta;
    // find beta based on the number of neighbor vertices
    if (n > 3) {
        beta = (3.0/(8.0 * n));
    }else if (n == 3) {
        beta = 3.0 / 16.0;
    }else if (n == 2) {
        beta = 1.0 / 8.0;
    }

    QVector3D new_point = QVector3D(0,0,0);
    // accumulated the neighbor results to get the new point
    for(long i = 0; i < (long) vertices_list[index].neighbor.size(); i++) {
        new_point += buffer[vertices_list[index].neighbor[i]] * beta;
    }

    // itself
    new_point += buffer[index] * (1-(beta * n));
    // update
    curr[index] = new_point;
}

void Mesh::loopSubdivision() {
    vector<QVector3D> vertices_buffer;
    vector<Mesh_Face> faces_buffer;

    // Firstly copy over everything in the vertices to the buffer
    for(uint i = 0; i < vertices.size(); i++) {
        vertices_buffer.push_back(vertices[i]);
    }

    // Secondly copy over all the faces to the buffer
    for(long i = 0; i < (long) faces.size(); i++) {
        faces_buffer.push_back(faces[i]);
    }

    long size = faces.size();
    for (long i = 0; i < (long) size; i++) {
        // for each face, we need to have three odd vertices
        // Each edge has one odd vertex
        // Iterate the neighbor faces

        // temporary array storing the original vertices
        vector<long> one_tri;
        one_tri.push_back(faces[i].vert[0]);
        one_tri.push_back(faces[i].vert[1]);
        one_tri.push_back(faces[i].vert[2]);

        for(long j = 0; j < 3; j++) {
            QVector3D mid = QVector3D(0,0,0);

            if (faces[i].adj[j] != -1) {
                // compute odd vertices
                // look for the edge which connect two vertices
                // look for the common and uncommon vertices from the current triangle
                for(long a = 0; a < (long) 3; a++) {
                    int count = 0;

                    for(long b = 0; b < (long) 3; b++) {
                        // common vertices
                        if (faces[i].vert[a] == faces[faces[i].adj[j]].vert[b]) {
                            mid += vertices_buffer[faces[i].vert[a]] * (3.0/8.0);
                            count++;
                        }
                    }
                    if (count == 0) {
                        // one of the uncommon vertex
                        mid += vertices_buffer[faces[i].vert[a]] * (1.0/8.0);
                    }
                }

                // look for the other uncommon vertex from the other triangle
                for(long a = 0; a < (long) 3; a++) {
                    int count = 0;
                    for(long b = 0; b < (long) 3; b++) {
                        // add up the count if seeing the common vertex
                        if (faces[faces[i].adj[j]].vert[a] == faces[i].vert[b]) {
                            count++;
                        }
                    }
                    if (count == 0) {
                        // the other uncommon vertex
                        mid += vertices_buffer[faces[faces[i].adj[j]].vert[a]] * (1.0/8.0);
                    }
                }

                // check if the vertex exists in the vertices
                int same = -1;
                for(uint n = 0; n < vertices.size(); n++) {
                    if (vertices[n] == mid) {
                        same = n;
                    }
                }

                if (same == -1) {
                    vertices.push_back(mid);
                    one_tri.push_back(vertices.size()-1);
                }else {
                    one_tri.push_back(same);
                }

            }else {
                // When the face has no adjacency face, just find the mid point as the new point
                long pos;
                if(j == 0) {
                    pos = insert_midPT(faces[i].vert[0],faces[i].vert[1],
                            vertices_buffer, vertices);
                }else if(j == 1) {
                    pos = insert_midPT(faces[i].vert[1],faces[i].vert[2],
                            vertices_buffer, vertices);
                }else if(j == 2) {
                    pos = insert_midPT(faces[i].vert[0],faces[i].vert[2],
                            vertices_buffer, vertices);
                }
                one_tri.push_back(pos);
            }
        }

        // compute even
        // calculate n and beta
        calculate_evenBeta(faces[i].vert[0], vertices_buffer, vertices);
        calculate_evenBeta(faces[i].vert[1], vertices_buffer, vertices);
        calculate_evenBeta(faces[i].vert[2], vertices_buffer, vertices);

        // make new faces with odd vertices
        Mesh_Face f1 = Mesh_Face(one_tri[0], one_tri[3], one_tri[5]);
        Mesh_Face f2 = Mesh_Face(one_tri[3], one_tri[1], one_tri[4]);
        Mesh_Face f3 = Mesh_Face(one_tri[2], one_tri[4], one_tri[5]);

        // updating the faces data structure
        faces_buffer[i].vert[0] = one_tri[3];
        faces_buffer[i].vert[1] = one_tri[4];
        faces_buffer[i].vert[2] = one_tri[5];
        faces_buffer.push_back(f1);
        faces_buffer.push_back(f2);
        faces_buffer.push_back(f3);
    }

    // rebuild mesh/ connect vertices to create new faces
    faces = faces_buffer;
}


void Mesh::add_face(const vector<int> &cur_vert) {
    if(cur_vert.size() > 3) {
        // If number of edges in face is greater than 3,
        // decompose into triangles as a triangle fan.
        int v0 = cur_vert[0], v1 = cur_vert[1], v2 = cur_vert[2];

        faces.push_back(Mesh_Face(v0, v1, v2));     // First face
        
        // all subsequent faces
        for( size_t i = 3; i < cur_vert.size(); i++ ) {
            v1 = v2; v2 = cur_vert[i];
            faces.push_back(Mesh_Face(v0, v1, v2));
        }
    }
    else if(cur_vert.size() == 3) {
        faces.push_back(Mesh_Face(cur_vert[0], cur_vert[1], cur_vert[2]));
    }
}

bool Mesh::load_obj(QString filename) {
    QFile objfile(filename);
    if (!objfile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return false; //error
    }
    QTextStream in(&objfile);
    long face_cnt = 0;

    while (!in.atEnd()) {
        QString line = in.readLine();
        line = line.trimmed();
        line = line.replace("\t", " ");

        QStringList tokens = line.trimmed().split(' ', QString::SkipEmptyParts);
        if(tokens.size() == 0) continue;

        if(tokens[0] == "v") {
            if(tokens.size() < 4) return false; // eror
            float x = tokens[1].toFloat();
            float y = tokens[2].toFloat();
            float z = tokens[3].toFloat();
            vertices.push_back(QVector3D(x,y,z));
        }
        if(tokens[0] == "f") {
            vector<int> cur_vert;
            for(int i = 1; i < tokens.size(); i++) {
                QStringList indexes = tokens[i].split("/");
                if(indexes.size() >= 1) {
                    if(indexes[0].toLong() < 0) {  cur_vert.push_back(vertices.size() + indexes[0].toLong()); }
                    else { cur_vert.push_back(indexes[0].toLong() - 1); }
                }
            }
            face_cnt++;
            add_face(cur_vert);
        }
    }
    cout << "face_cnt=" << face_cnt << endl;
    cout << "faces.size()=" << faces.size() << endl;
    cout << "vertices.size()=" << vertices.size() << endl;

    recenter();
    return true;
}

void Mesh::recenter() {
    if( vertices.size() < 1) return;
    QVector3D maxPoint = vertices[0];
    QVector3D minPoint = vertices[0];

    // Find the AABB
    for( uint i = 0; i < vertices.size(); ++i ) {
        QVector3D & point = vertices[i];
        if( point[0] > maxPoint[0] ) maxPoint[0] = point[0];
        if( point[1] > maxPoint[1] ) maxPoint[1] = point[1];
        if( point[2] > maxPoint[2] ) maxPoint[2] = point[2];
        if( point[0] < minPoint[0] ) minPoint[0] = point[0];
        if( point[1] < minPoint[1] ) minPoint[1] = point[1];
        if( point[2] < minPoint[2] ) minPoint[2] = point[2];
    }

    // Center of the AABB
    QVector3D center = QVector3D( (maxPoint[0] + minPoint[0]) / 2.0f,
    (maxPoint[1] + minPoint[1]) / 2.0f,
    (maxPoint[2] + minPoint[2]) / 2.0f );

    // Translate center of the AABB to the origin
    for( uint i = 0; i < vertices.size(); ++i ) {
        QVector3D & point = vertices[i];
        point = point - center;
    }
}


void Mesh::process_example() {
    for(size_t v = 0; v < vertices.size(); v++) {
        if(vertices[v][0] > 0) {
            vertices[v][0] += 3.5;
        }
    }
}

void Mesh::storeVBO() {
    vector<QVector3D> tri_vert, tri_bary;

    for(long f = 0; f < (long)faces.size(); f++) {
        tri_vert.push_back(vertices.at(faces[f].vert[0]));
        tri_vert.push_back(vertices.at(faces[f].vert[1]));
        tri_vert.push_back(vertices.at(faces[f].vert[2]));

        tri_bary.push_back(QVector3D(1,0,0));;
        tri_bary.push_back(QVector3D(0,1,0));;
        tri_bary.push_back(QVector3D(0,0,1));;
    }

    if(vertexBuffer.isCreated()) vertexBuffer.destroy();
    vertexBuffer.create();
    vertexBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
    vertexBuffer.bind();
    vertexBuffer.allocate(&tri_vert[0] , sizeof( QVector3D ) * tri_vert.size());

    if(baryBuffer.isCreated()) baryBuffer.destroy();
    baryBuffer.create();
    baryBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
    baryBuffer.bind();
    baryBuffer.allocate(&tri_bary[0] , sizeof( QVector3D ) * tri_bary.size());
}
