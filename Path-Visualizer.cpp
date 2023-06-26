// VTK
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTubeFilter.h>
#include <vtkTransform.h>
#include <vtkTextProperty.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkUnsignedCharArray.h>
#include <vtkNamedColors.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkCellData.h>
#include <vtkActorCollection.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>


// STD
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <array>
using namespace std;

int nbr_paths{ 0 };
int nbr_actors{ 0 };

bool visibility = true;
vector<bool> vis;
bool centerline_visibility = true;
bool normal_visibility = true;
bool solid_visibility = true;
bool clip = false;

// Colors
double lime_green[3]{ 50, 205, 50 };  // Perimeter
double cyan[3]{ 0, 255, 255 }; // Skeleton
double red[3]{ 255, 0, 0 };    // Arrows
map<int, double(*)[3]> color_map{ pair<int, double(*)[3]>(1, &lime_green), pair<int, double(*)[3]>(16, &cyan) };

// Define interaction style
class InteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static InteractorStyle* New();
    vtkTypeMacro(InteractorStyle, vtkInteractorStyleTrackballCamera);

    // Display the coordinates of the closest point
    // Should only be used with centerlines
    virtual void OnRightButtonDown() override
    {
        this->Interactor->GetPicker()->Pick(
            this->Interactor->GetEventPosition()[0],
            this->Interactor->GetEventPosition()[1],
            0, // always zero.
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

        double picked[3];
        this->Interactor->GetPicker()->GetPickPosition(picked);
        std::cout << fixed << setprecision(6) << "Point(" << picked[0] << ", " << picked[1] << ", " << picked[2] << ")" << std::endl;

        // Forward events
        vtkInteractorStyleTrackballCamera::OnRightButtonDown();
    }

    // Keyboard
    virtual void OnKeyPress() override
    {
        // Get the keypress
        vtkRenderWindowInteractor* rwi = this->Interactor;
        std::string key = rwi->GetKeySym();

        // Rotate CW
        if (key == "Left")
        {
            this->CurrentRenderer->GetActiveCamera()->SetViewUp(0, 0, 1);
            this->CurrentRenderer->GetActiveCamera()->Azimuth(-1);
            this->CurrentRenderer->GetRenderWindow()->Render();
        }

        // Rotate CCW
        if (key == "Right")
        {
            this->CurrentRenderer->GetActiveCamera()->SetViewUp(0, 0, 1);
            this->CurrentRenderer->GetActiveCamera()->Azimuth(1);
            this->CurrentRenderer->GetRenderWindow()->Render();
        }

        // Toggle visibility of all paths
        if (key == "v")
        {
            visibility = !visibility;
            centerline_visibility = visibility;
            normal_visibility = visibility;
            solid_visibility = visibility;

            for (int i = 0; i < nbr_actors; ++i)
            {
                auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(i));
                act->SetVisibility(visibility);
                this->Interactor->Render();
            }

            std::fill(vis.begin(), vis.end(), visibility);
        }

        // Toggle visible paths' centerline visibility
        if (key == "c")
        {
            centerline_visibility = !centerline_visibility;

            for (int i = 0; i < nbr_paths; ++i)
            {
                if (vis[i])
                {
                    auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(i));
                    act->SetVisibility(centerline_visibility);
                    this->Interactor->Render();
                }
            }
        }

        // Toggle visible paths' normal visibility
        if (key == "n")
        {
            normal_visibility = !normal_visibility;

            for (int i = 0; i < nbr_paths; ++i)
            {
                if (vis[i])
                {
                    auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(nbr_paths + i));
                    act->SetVisibility(normal_visibility);
                    this->Interactor->Render();
                }
            }
        }

        // Toggle visible paths' width visibility
        if (key == "s")
        {
            solid_visibility = !solid_visibility;

            for (int i = 0; i < nbr_paths; ++i)
            {
                if (vis[i])
                {
                    auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(nbr_paths * 2 + i));
                    act->SetVisibility(solid_visibility);
                    this->Interactor->Render();
                }
            }
        }

        // Make individual paths visible/invisible
        if (std::isdigit(key.front()))
        {
            int k = std::stoi(key);

            if (k < nbr_paths)
            {
                vis[k] = !vis[k];

                if (!centerline_visibility && !normal_visibility && !solid_visibility)
                    centerline_visibility = true;

                for (int i = k; i < nbr_actors; i += nbr_paths)
                {
                    if (centerline_visibility && i < nbr_paths)
                    {
                        auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(i));
                        act->SetVisibility(vis[k]);
                        this->Interactor->Render();
                    }

                    if (normal_visibility && i >= nbr_paths && i < 2 * nbr_paths)
                    {
                        auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(i));
                        act->SetVisibility(vis[k]);
                        this->Interactor->Render();
                    }

                    if (solid_visibility && i >= 2 * nbr_paths)
                    {
                        auto act = vtkActor::SafeDownCast(this->CurrentRenderer->GetActors()->GetItemAsObject(i));
                        act->SetVisibility(vis[k]);
                        this->Interactor->Render();
                    }
                }
            }
        }
    }
};
vtkStandardNewMacro(InteractorStyle);

string dir = "Pathing_Visualization_Files";
void import_data(vector<vector<array<double, 3>>>& points, vector<vector<array<double, 3>>>& normals, vector<vector<double>>& widths, vector<vector<int>>& types)
{
    // Points
    ifstream pfs(dir + "/points.dat", ifstream::binary);
    istringstream pss(ios_base::binary);
    string point_line;
    while (getline(pfs, point_line, '\n'))
    {
        vector<array<double, 3>> path_points;
        array<double, 3> point;

        pss.str(point_line);
        while (pss >> point[0] && pss >> point[1] && pss >> point[2])
            path_points.push_back(point);
        pss.clear();

        points.push_back(path_points);
    }
    pfs.close();


    // Normals
    ifstream nfs(dir + "/normals.dat", ifstream::binary);
    istringstream nss(ios_base::binary);
    string normal_line;
    while (getline(nfs, normal_line, '\n'))
    {
        vector<array<double, 3>> path_normals;
        array<double, 3> normal;

        nss.str(normal_line);
        while (nss >> normal[0] && nss >> normal[1] && nss >> normal[2])
            path_normals.push_back(normal);
        nss.clear();

        normals.push_back(path_normals);
    }
    nfs.close();


    // Widths
    ifstream wfs(dir + "/widths.dat", ifstream::binary);
    istringstream wss(ios_base::binary);
    string width_line;
    while (getline(wfs, width_line, '\n'))
    {
        double w;
        vector<double> path_widths;

        wss.str(width_line);
        while (wss >> w)
            path_widths.push_back(w);
        wss.clear();

        widths.push_back(path_widths);
    }
    wfs.close();


    // Types
    ifstream tfs(dir + "/types.dat", ifstream::binary);
    istringstream tss(ios_base::binary);
    string type_line;
    while (std::getline(tfs, type_line, '\n'))
    {
        int t;
        vector<int> path_types;

        tss.str(type_line);
        while (tss >> t)
            path_types.push_back(t);
        tss.clear();

        types.push_back(path_types);
    }
    tfs.close();
}

int main(int, char* [])
{
    vtkObject::GlobalWarningDisplayOff();

    // Import data
    vector<vector<array<double, 3>>> points, normals;
    vector<vector<double>> widths;
    vector<vector<int>> types;
    import_data(points, normals, widths, types);

    nbr_paths = points.size();
    nbr_actors = nbr_paths * 3;
    vis = vector<bool>(nbr_paths, true);

    vector<vtkSmartPointer<vtkPolyData>> centerlines;
    vector<vtkSmartPointer<vtkTubeFilter>> tubes;
    vector<vtkSmartPointer<vtkGlyph3D>> arrows;

    for (int path_idx = 0; path_idx < nbr_paths; ++path_idx)
    {
        vtkNew<vtkPoints> vtk_points;

        vtkNew<vtkCellArray> vtk_cells;
        vtk_cells->InsertNextCell(points[path_idx].size());

        vtkNew<vtkDoubleArray> vtk_normals;
        vtk_normals->SetName("normals");
        vtk_normals->SetNumberOfComponents(3);

        vtkNew<vtkDoubleArray> vtk_radii;
        vtk_radii->SetName("radii");

        vtkNew<vtkUnsignedCharArray> path_colors;
        path_colors->SetName("Path_Colors");
        path_colors->SetNumberOfComponents(3);

        vtkNew<vtkUnsignedCharArray> arrow_colors;
        arrow_colors->SetName("Arrow_Colors");
        arrow_colors->SetNumberOfComponents(3);

        // Construct path
        for (int segment_idx = 0; segment_idx < types[path_idx].size(); ++segment_idx)
        {
            // Insert points
            auto id = vtk_points->InsertNextPoint(points[path_idx][segment_idx].data());
            vtk_cells->InsertCellPoint(id);

            // Insert normals
            vtk_normals->InsertNextTuple(normals[path_idx][segment_idx].data());

            // Insert radii
            vtk_radii->InsertNextValue(widths[path_idx][segment_idx] / 2.0);

            // Insert path color
            path_colors->InsertNextTuple(*color_map[types[path_idx][segment_idx]]);

            // Insert arrow color
            arrow_colors->InsertNextTuple(red);
        }

        // Construct end of path
        auto id = vtk_points->InsertNextPoint(points[path_idx].back().data());
        vtk_cells->InsertCellPoint(id);
        vtk_radii->InsertNextValue(widths[path_idx].back() / 2.0);
        vtk_normals->InsertNextTuple(normals[path_idx].back().data());
        path_colors->InsertNextTuple(*color_map[types[path_idx].back()]);
        arrow_colors->InsertNextTuple(red);

        // Construct polyData to hold path
        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(vtk_points);
        polyData->SetLines(vtk_cells);
        polyData->GetPointData()->AddArray(vtk_normals);
        polyData->GetPointData()->SetActiveVectors("normals");
        polyData->GetPointData()->AddArray(vtk_radii);
        polyData->GetPointData()->SetActiveScalars("radii");
        polyData->GetPointData()->AddArray(path_colors);
        polyData->GetPointData()->AddArray(arrow_colors);


        if (clip)
        {
            double center[3];
            if (centerlines.empty())
                polyData->GetCenter(center);

            vtkNew<vtkPlane> clip_plane;
            clip_plane->SetOrigin(center);
            clip_plane->SetNormal(-1.0, 0.0, 0.0);
            vtkNew<vtkClipPolyData> clipper;
            clipper->SetClipFunction(clip_plane);
            clipper->SetInputData(polyData);
            clipper->GenerateClippedOutputOn();
            clipper->Update();
            polyData = clipper->GetClippedOutput();
        }

        // Create centerlines of paths
        centerlines.push_back(polyData);

        // Create tube filters for polyData
        vtkNew<vtkTubeFilter> tube;
        tube->SetInputData(polyData);
        tube->SetNumberOfSides(20);
        tube->CappingOn();
        tube->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
        tubes.push_back(tube);

        // Create arrows for polyData
        vtkNew<vtkArrowSource> arrow;
        vtkNew<vtkGlyph3D> glyph;
        glyph->SetSourceConnection(arrow->GetOutputPort());
        glyph->SetInputData(polyData);
        glyph->SetColorMode(1);
        glyph->SetScaleFactor(1);
        glyph->SetScaleModeToScaleByVector();
        glyph->SetColorModeToColorByScalar();
        glyph->Update();
        arrows.push_back(glyph);
    }


    // Set up centerline mapper and actors
    vector<vtkSmartPointer<vtkActor>> centerline_actors;
    for (vtkSmartPointer<vtkPolyData> centerline : centerlines)
    {
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(centerline);
        mapper->ScalarVisibilityOn();
        mapper->SetScalarModeToUsePointFieldData();
        mapper->SelectColorArray("Path_Colors");

        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        centerline_actors.push_back(actor);
    }

    // Set up tube mapper and actors
    vector<vtkSmartPointer<vtkActor>> tube_actors;
    for (vtkSmartPointer<vtkTubeFilter> tube : tubes)
    {
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(tube->GetOutputPort());
        mapper->ScalarVisibilityOn();
        mapper->SetScalarModeToUsePointFieldData();
        mapper->SelectColorArray("Path_Colors");

        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        tube_actors.push_back(actor);
    }

    // Set up arrow mapper and actors
    vector<vtkSmartPointer<vtkActor>> arrow_actors;
    for (vtkSmartPointer<vtkGlyph3D> arrow : arrows)
    {
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(arrow->GetOutputPort());
        mapper->SetScalarModeToUsePointFieldData();
        mapper->SelectColorArray("Arrow_Colors");

        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        arrow_actors.push_back(actor);
    }


    // Setup render window, renderer, and interactor
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    vtkNew<InteractorStyle> style;
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindow->SetWindowName("Geomerty Visualizer");
    renderWindow->AddRenderer(renderer);
    renderWindowInteractor->SetInteractorStyle(style);
    vtkNew<vtkPointPicker> pointPicker;
    renderWindowInteractor->SetPicker(pointPicker);
    style->SetCurrentRenderer(renderer);
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderer->SetBackground(vtkNew<vtkNamedColors>()->GetColor3d("DarkSlateGray").GetData());


    // Add centerline actors to renderer
    for (auto centerline_actor : centerline_actors)
        renderer->AddActor(centerline_actor);

    // Add arrow actors to renderer
    for (auto arrow_actor : arrow_actors)
        renderer->AddActor(arrow_actor);

    // Add tube actors to renderer
    for (auto tube_actor : tube_actors)
        renderer->AddActor(tube_actor);


    // Camera
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Elevation(-90); //(-60);
    renderer->GetActiveCamera()->SetViewUp(0, 0, 1);
    renderer->GetActiveCamera()->Azimuth(-46);
    double clip[2]{ 0.1, 1000.0 };
    renderer->GetActiveCamera()->SetClippingRange(clip);
    renderer->GetActiveCamera()->SetViewAngle(65);
    renderer->GetActiveCamera()->Dolly(2.9);
    
    
    

    // Render
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}