#ifndef __GLVIEW_HPP__
#define __GLVIEW_HPP__

#include <QtGui>
#include <QOpenGLWidget>
#include "Mesh.hpp"
#include <iostream>

class GLview : public QOpenGLWidget, protected QOpenGLFunctions  {
    Q_OBJECT
public:
    // constructor. DO NOT MODIFY
    GLview(QWidget *parent = 0);

    // destructor. DO NOT MODIFY
    ~GLview();

    void keyPressGL(QKeyEvent* e);

    bool LoadOBJFile(const QString file);

protected:
    void initializeGL(); // QGLWidget OpenGL interface
    void initCameraGL();
    void paintGL();
    void resizeGL(int width, int height);

    Mesh *mesh;  // Geometry data.
    QOpenGLVertexArrayObject vao; // vertex buffer
    QOpenGLShaderProgram wire_shader;

    // Compiles the shader programs. DO NOT MODIFY.
    bool prepareShaderProgram(QOpenGLShaderProgram &shader, const QString &vertex_file, const QString &fragment_file);
    bool prepareShaderProgramTex();

    void update_mesh();

    // Camera parameters --------------
    QVector3D eye, lookCenter, lookUp;
    QQuaternion camrot;
    float yfov, neardist, fardist;
    // -------------------------------

    // Used to track mouse position for keyboard shortcuts
    // prior to mouse driven interaction.
    bool lastPosFlag; // Flag is true if a last position has been captured. Set to false
    float lastPosX, lastPosY;   // Mouse state information.
    bool scaleFlag, translateFlag, rotateFlag;  // Camera movement state information.

    // Store mouse press position. DO NOT MODIFY
    void mousePressEvent(QMouseEvent *event);

    // Update camera position on mouse movement. DO NOT MODIFY.
    void mouseMoveEvent(QMouseEvent *event);

    // DO NO MODIFY
    void toggleRotate();

    // DO NO MODIFY
    void toggleScale();

    // DO NO MODIFY
    void toggleTranslate();
private:
    void ShowContextMenu(const QMouseEvent *event);

private slots:
    // implement these operations
    void process_example(); // this is an example
    void inflate();
    void randomNoise();
    void smooth();
    void sharpen();
    void truncate();
    void splitFaces();
    void starFaces();
    void splitLongEdges();
    void collapseShortEdges();
    void loopSubdivision();
    void centerVerticesTangentially();
    void flipEdges();
    void crop();
    void bilateralSmoothing();
    void meshSimplification();

};

#endif
