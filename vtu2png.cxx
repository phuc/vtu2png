#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPNGWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkRenderer.h>
#include <vtkGraphicsFactory.h>
#include <vtkImagingFactory.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAppendPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkScalarsToColors.h>
#include <vtkLookupTable.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkLegendBoxActor.h>
#include <vtkSphereSource.h>
int main (int argc, char *argv[])
{
  if (argc < 3)
    {
    std::cerr << "Usage: " << argv[0]
              << " Input(.vtu) Magnification Scalar Output(.png) [Magnification]"
              << std::endl;
    return EXIT_FAILURE;
    }

  int magnification = 4;
  // Setup offscreen rendering
  vtkSmartPointer<vtkGraphicsFactory> graphics_factory = 
    vtkSmartPointer<vtkGraphicsFactory>::New();
  graphics_factory->SetOffScreenOnlyMode( 1);
  graphics_factory->SetUseMesaClasses( 1 );

  vtkSmartPointer<vtkImagingFactory> imaging_factory = 
    vtkSmartPointer<vtkImagingFactory>::New();
  imaging_factory->SetUseMesaClasses( 1 );
  
  // unstructured grid reader    
  vtkSmartPointer<vtkXMLUnstructuredGridReader> unstructuredReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  unstructuredReader->SetFileName(argv[1]);
  unstructuredReader->Update();
      
  // convert
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();  
  surfaceFilter->SetInputConnection(unstructuredReader->GetOutputPort());
  surfaceFilter->Update();  
      
  double scalarRange[2];
  
  vtkSmartPointer<vtkPointData> pointData =  vtkSmartPointer<vtkPointData>::New();
  pointData = surfaceFilter->GetOutput()->GetPointData();
  pointData->SetActiveScalars(pointData->GetArrayName(0));
  pointData->GetScalars()->GetRange(scalarRange);

  vtkSmartPointer<vtkAppendPolyData> appendFilledContours = vtkSmartPointer<vtkAppendPolyData>::New();

  int numberOfContours = 24;

  double delta = (scalarRange[1] - scalarRange[0]) / static_cast<double> (numberOfContours - 1);
  
  // Keep the clippers alive
  std::vector<vtkSmartPointer<vtkClipPolyData> > clippersLo;
  std::vector<vtkSmartPointer<vtkClipPolyData> > clippersHi;

  for (int i = 0; i < numberOfContours; i++)
    {
    double valueLo = scalarRange[0] + static_cast<double> (i) * delta;
    double valueHi = scalarRange[0] + static_cast<double> (i + 1) * delta;
    clippersLo.push_back(vtkSmartPointer<vtkClipPolyData>::New());
    clippersLo[i]->SetValue(valueLo);
    if (i == 0)
      {
      clippersLo[i]->SetInputConnection(surfaceFilter->GetOutputPort());
      }
    else
      {
      clippersLo[i]->SetInputConnection(clippersHi[i - 1]->GetOutputPort(1));
      }
    clippersLo[i]->InsideOutOff();
    clippersLo[i]->Update();

    clippersHi.push_back(vtkSmartPointer<vtkClipPolyData>::New());
    clippersHi[i]->SetValue(valueHi);
    clippersHi[i]->SetInputConnection(clippersLo[i]->GetOutputPort());
    clippersHi[i]->GenerateClippedOutputOn();
    clippersHi[i]->InsideOutOn();
    clippersHi[i]->Update();
    if (clippersHi[i]->GetOutput()->GetNumberOfCells() == 0)
      {
      continue;
      }
    
    vtkSmartPointer<vtkFloatArray> cd =
      vtkSmartPointer<vtkFloatArray>::New();
    cd->SetNumberOfComponents(1);
    cd->SetNumberOfTuples(clippersHi[i]->GetOutput()->GetNumberOfCells());
    cd->FillComponent(0, valueLo);

    clippersHi[i]->GetOutput()->GetCellData()->SetScalars(cd);
    appendFilledContours->AddInputConnection(clippersHi[i]->GetOutputPort());
    }

  vtkSmartPointer<vtkCleanPolyData> filledContours =
    vtkSmartPointer<vtkCleanPolyData>::New();
  filledContours->SetInputConnection(appendFilledContours->GetOutputPort());

  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetTableRange(scalarRange[0], scalarRange[1]);
  lut->SetNumberOfTableValues(numberOfContours + 1);
  lut->SetHueRange(0.667, 0.0);
  lut->Build();
  
  vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  contourMapper->SetInputConnection(filledContours->GetOutputPort());
  contourMapper->SetScalarRange(scalarRange[0], scalarRange[1]);
  contourMapper->SetScalarModeToUseCellData();
  contourMapper->ScalarVisibilityOn();
  contourMapper->SetColorModeToMapScalars();
  contourMapper->SetColorMode(2);
  contourMapper->SetLookupTable(lut);
 
  vtkSmartPointer<vtkActor> contourActor =
    vtkSmartPointer<vtkActor>::New();
  contourActor->SetMapper(contourMapper);
  contourActor->GetProperty()->SetInterpolationToFlat();

  vtkSmartPointer<vtkContourFilter> contours =
    vtkSmartPointer<vtkContourFilter>::New();
  contours->SetInputConnection(filledContours->GetOutputPort());
  contours->GenerateValues(numberOfContours, scalarRange[0], scalarRange[1]);
  
  vtkSmartPointer<vtkLegendBoxActor> legend = vtkSmartPointer<vtkLegendBoxActor>::New();
  legend->SetNumberOfEntries(contours->GetNumberOfContours());
  
  for (int i = 0; i < contours->GetNumberOfContours(); i++) {
    char text[128];
    sprintf(text, "%.2f", contours->GetValue(contours->GetNumberOfContours() - 1 - i));    
    vtkSmartPointer<vtkSphereSource> legendSphereSource = vtkSmartPointer<vtkSphereSource>::New();
    vtkSmartPointer<vtkPolyData> legendSphere = legendSphereSource->GetOutput();
    legend->SetEntry(i, legendSphere, text, lut->GetTableValue(contours->GetNumberOfContours() - 1 - i));
  }
    
  /*
  // The usual renderer, render window and interactor
  vtkSmartPointer<vtkRenderer> ren1 =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor>
    iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  ren1->SetBackground(.1, .2, .3);
  renWin->AddRenderer(ren1);
  iren->SetRenderWindow(renWin);

  // Add the actors
  ren1->AddActor(contourActor);
  ren1->AddActor(contourLineActor);

  // Begin interaction
  renWin->Render();
  iren->Start();
  */
  
  //======= Offscreen Render 
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(contourActor);
  renderer->AddActor(legend);
  
  // render Window + renderer
  vtkSmartPointer<vtkRenderWindow> renderWindow = 
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetOffScreenRendering( 1 ); 
  renderWindow->AddRenderer(renderer);

  // Let the renderer compute good position and focal point
  //renderer->GetActiveCamera()->Azimuth(30);
  //renderer->GetActiveCamera()->Elevation(30);
  //renderer->ResetCamera();
  //renderer->GetActiveCamera()->Dolly(1.4);
  //renderer->ResetCameraClippingRange();
  renderer->SetBackground(.3, .4, .5);

  renderWindow->SetSize(1024, 1024);  
  renderWindow->Render();
  
  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = 
      vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetMagnification(magnification);
  windowToImageFilter->Update();    
  
  //vtkSmartPointer<vtkRenderLargeImage> renderLarge =
  //  vtkSmartPointer<vtkRenderLargeImage>::New();
  //renderLarge->SetInput(renderer);
  //renderLarge->SetMagnification(magnification);

  std::cout << "Saving image in " << argv[2] << std::endl;
  vtkSmartPointer<vtkPNGWriter> writer =
    vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(argv[2]);
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
  
  return EXIT_SUCCESS;
}
