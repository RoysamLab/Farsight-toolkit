/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

/*=========================================================================
  Program:   Open Snake Tracing System
  Autohr:    Yu Wang
  Email: wangy15@rpi.edu
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $
=========================================================================*/

#ifndef TRACINGGUI_H
#define TRACINGGUI_H

#include <QtGui>
#include "TracingView.h"
#include "ParametersGroup.h"
#include "SlidersGroup.h"
#include "TracingCore/Montage/Register.h"
#include "dialogs.h"
#include "dialogs_montage.h"
#include "time.h"
#include <vnl/vnl_random.h> 
#include "Rendering/ImageActors.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"


QImage convertITK2QT(ImagePointer2D I, bool mono);
QImage convertITK2QT(ProbImagePointer2D I, bool mono);
QImage convertITK2QT(LabelImagePointer2D I, bool mono);
QImage convertITK2QT(RGBImagePointer2D I);


class QtTracer : public QMainWindow
{
Q_OBJECT

public:
	QtTracer(QWidget * parent = 0, Qt::WindowFlags flags = 0);

public slots:
	void exit();
	void about();
	void stackImage();
	void loadImage();
	void loadDisplayImage();
	void reloadImage();

	void autoProcessing();

	void batchPreprocessing();

	void batchProcessingI();
	void batchTracingI();

	void batchProcessingII();
	void generateGraph();
	void createMontage();
	void batchTracingII();

	void addGlobalSnake();
	int estimateZShift();

	void outputSeeds();

	void set_cross_section_view(int in);
	void tab_change(int in);

	void tracing_suspend();
	void tracing_go();
    void tracing_stop();

	void outsidePreprocess();
    void curvelet_scalar_voting();

	void Preprocess();
	void resetPreprocessStage();
	void Reprocess_Vessel(int in);

	void Init_Tracing();
	void Tracing();
	void Seed_Snake_Tracing();
    void Manual_Seed_Tracing(PointList3D seeds);
    void Manual_Seed_TracingI(PointList3D seeds);
    void SegmentationI();
	void Radius_Estimation();
	void removeIsolated();
	void enable_manual_seeding();

	void Process();
	void nextStep();
	void backStep();
	void changeDisplay(int in);
	void changeSlice(int in);
	void changeOpacity(int in);
	void changeOpacityTh(int in);

    void SnakesChangedSlot();
	void PointClickedSlot(Point3D pt);
   
	void output_swc_raw_slot();
	void outputSWC();
    void outputSWC_Raw();
	void outputSWC_Raw_MT();
    void outputVTK();

    void setRootPoint();
	void convertSnakeTree();
    void clearSnakeTree();
	void findBranch(int snake_label, int root_label, Point3D root_point, vnl_vector<int>*snake_visit_label, int *point_id, SnakeTree *tree);
    void findBranch_Raw(int snake_label, int root_label, Point3D root_point, vnl_vector<int>*snake_visit_label, int *point_id, QTextStream *out_txt, PointList3D *wrote_pt, PointList3D *all_pt, std::vector<int> *bl);
	void findBranch_Raw_MT(int snake_label, int root_label, Point3D root_point, vnl_vector<int>*snake_visit_label, int *point_id, QTextStream *out_txt, PointList3D *wrote_pt, PointList3D *all_pt, std::vector<int> *bl);
	void refineBranch();
	void pickSomaSeeds();
	void segmentSoma();
	void clearSegmentation();
	void saveSoma();
	void removeSoma();
	void saveImage();

	void save_settings();
	void load_settings();

	void auto_mode_slot(bool in);
	void batchI_mode_slot(bool in);
	void batchII_mode_slot(bool in);

	void draw3DTraces();
	void draw3DTraces_Global();
	void drawVolume();
	void draw3DTracing(SnakeClass s);
	void surfaceRendering(bool rand_color);
	void surfaceRendering(bool rand_color, ImagePointer ID);
	void writeRendering(vtkRenderWindow *rw, const char *filename);

    //edits
    void deleteSnake();
	void mergeSnake();
	void splitSnake();
	void branchSnake();
	void updateFN(float fn);
	void breakBranches();

	//slots for vtk events
	void vtk_right_pick(vtkObject * obj);
	void vtk_left_pick(vtkObject * obj);
	void vtk_removePoint();

signals:
	void selected(const QItemSelection &, QItemSelectionModel::SelectionFlags command);
	void sendSelected(int * sr);

protected:
    void closeEvent();

private:
    
	//Widgets
    QScrollArea *scrollArea, *scrollArea_BW, *scrollArea_MT;
	QSplitter *main_splitter;
    QDockWidget *dockWidget;
	QWidget *right_widget, *left_widget;
    SlidersGroup *slider;
	QPushButton *backButton, *nextButton, *startButtonI, *startButtonII, *quitButton;

    QTabWidget *imageTab;
	QLabel *imageLabel;
	QLabel *crossSectionLabel1;
	QLabel *crossSectionLabel2;
	QWidget *rendering_widget;

	QLabel *montageLabel;

	//labels showing the length precision and recall rate, and number of edits
	QLabel *PrecisionLabel, *RecallLabel, *BranchLabel, *SplitLabel, *DeleteLabel, *MergeLabel;
    EditValidation *edits;

	QVBoxLayout *layout1, *layout2;
	QHBoxLayout *layout_paras, *layout3;
	QGridLayout *grid_layout;
    
    void createMenus();
	void createActions();
	void createMainGUI();
	void createToolBar();
	void initializeRendering();

	//Custom Widgets
	TracingView *tracingViewer;
    
	QToolBar *Tool;
	QToolBar *Tool_Misc;
	QMenu *fileMenu;
	QMenu *SeedMenu;
	QMenu *SnakeMenu;
	QMenu *BranchMenu;
	QMenu *TreeMenu;
	QMenu *SomaMenu;
	QMenu *PostSegMenu;
	QMenu *OutlierMenu;
 	QMenu *viewMenu;
	QMenu *aboutMenu;
    
	QAction *zoomInAct;
	QAction *zoomOutAct;
	QAction *normalSizeAct;
	QAction *aboutAct;
	QAction *loadImageAct;
	QAction *exitAct;

	QProgressBar *Progress;
	QAction *Go;
	QAction *Pause;
    QAction *Stop;
	QAction *Pick_Seed;
	QAction *Delete_Seed;
	QAction *Invert_Selection;
    QAction *Delete_Snake;
	QAction *Split_Snake;
	QAction *Merge_Snake;
	QAction *Set_Root;
	QAction *Clear_Tree;
	QAction *Create_Branch;
    QAction *Refine_Branch;
	QAction *Remove_Soma;
	QAction *Pick_Soma_Seeds;
	QAction *Segment_Soma;
	QAction *Save_Soma;
	QAction *Clear_Segmentation;
	QAction *SegmentI;
	QAction *SegmentII;
	QAction *Remove_Isolated;
	QAction *Estimate_Radius;
	QAction *Save_Image;
	QAction *Save_Setting;
	QAction *Load_Setting;
	QAction *Break_Branches;

	ImageOperation *IM;
    OpenSnakeTracer *Tracer; 

	General_Parameters *general_para;
	General_Parameters1 *general_para1; 
	General_Parameters12 *general_para12;
	General_Parameters2 *general_para2;
	General_Parameters3 *general_para3;
   
	int Process_Stage;
	int Old_Process_Stage;
	int Preprocess_Stage;

	bool tracing_started;

	time_t Start_t, End_t; 
	int preprocessing_time;
	int tracing_time;

	bool seeding;
  
	bool automated_process;
	bool batch_processI;
	bool batch_processII;
	bool batch_preprocessing;

	bool soma_tracing;

	bool using_stacked_image;

	bool outside_preprocess;
  
	int thread_finished;

	QString outside_preprocessor_path;

	PointList3D m_seed;

	SnakeTree *snake_tree;
	SnakeListClass seed_snakes;

    QString file, file_raw, file_display;

	setRootDialog *set_root_dialog;
	setTransformDialog *set_transform_dialog;

	QDir images_dir;
	QString swcs_path;
	int current_idx;

	int original_seed_num;

	QString *string_array; //the names of all images in the folder
    QString txt_fileName; //the txt file containing coordinates

	Register reg;
	PointList3D coordinates;
	bool montage_created;
	bool bpII_finished;
	SnakeListClass GSnakeList, GSnakeList_BP, temp_snake_list; //global snake list
	int current_lvl; //current level in the tile graph
	float min_x;
	int Z_shift; //estimated shift in Z dimension

	QImage Background; //montage image
	QDir xmls_dir;

    PointList3D wrote_pt, All_Pt;

	QSlider *opacity_slider;
    QLabel *opacity;
	QSlider *opacity_th_slider;
	QLabel *opacity_th;
	//Rendering
	QVTKWidget *QVTK;	
	ImageRenderActors *ImageActors;
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkCamera> Camera;

	//vtkSmartPointer<vtkActor> line_actor;
	vtkSmartPointer<vtkActor> point_actor;
	vtkSmartPointer<vtkActor> seed_actor;
	vtkSmartPointer<vtkBoundingBox> BoundingBox;

	std::vector< vtkSmartPointer<vtkActor> > line_actors;
	std::vector< vtkSmartPointer<vtkActor> > tube_actors;
    vtkSmartPointer<vtkActor> tracing_line_actor;

	std::vector< vtkSmartPointer<vtkActor> > branch_actors;
	std::vector< vtkSmartPointer<vtkActor> > boundary_actors;
	std::vector< vtkSmartPointer<vtkActor> > mesh_actors;
	std::vector< vtkSmartPointer<vtkActor> > soma_mesh_actors;

	vnl_matrix<double> soma_color;

	//vtk event to Qt Slot
	vtkEventQtSlotConnect *Connections;
	std::vector< vtkSmartPointer<vtkActor> > pick_sphere_actors;
	PointList3D picked_pts;
	bool soma_seeding;

};

#endif
