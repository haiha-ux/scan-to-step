import sys, numpy as np, trimesh, vtk
from pathlib import Path
out = Path('solution/img'); out.mkdir(exist_ok=True)
for i in (1,2,3):
    m=trimesh.load(f'scan_{i}.stl'); to,ext=trimesh.bounds.oriented_bounds(m); m.apply_transform(to)
    f=out/f'_obb_{i}.stl'; m.export(f)
    for nm,pos,up in (('iso',(1,-1.3,0.9),(0,0,1)),('iso2',(-1,1.2,-0.8),(0,0,1))):
        r=vtk.vtkSTLReader(); r.SetFileName(str(f)); r.Update()
        mp=vtk.vtkPolyDataMapper(); mp.SetInputConnection(r.GetOutputPort())
        a=vtk.vtkActor(); a.SetMapper(mp); a.GetProperty().SetColor(0.8,0.8,0.75)
        ren=vtk.vtkRenderer(); ren.SetBackground(0.3,0.33,0.38); ren.AddActor(a)
        c=ren.GetActiveCamera(); c.SetPosition(*[200*p for p in pos]); c.SetFocalPoint(0,0,0); c.SetViewUp(*up); ren.ResetCamera()
        w=vtk.vtkRenderWindow(); w.SetOffScreenRendering(1); w.AddRenderer(ren); w.SetSize(800,600); w.Render()
        g=vtk.vtkWindowToImageFilter(); g.SetInput(w); g.Update()
        p=vtk.vtkPNGWriter(); p.SetFileName(str(out/f'scan_{i}_{nm}.png')); p.SetInputConnection(g.GetOutputPort()); p.Write()
    print(i, np.round(ext,3))
