from node_runner.dialogs.editors import (
    MaterialEditorDialog, PropertyEditorDialog, ElementEditorDialog,
)
from node_runner.dialogs.creators import (
    CreateMaterialDialog, CreatePropertyDialog,
    CreateNodesDialog, CreateLineElementDialog, CreatePlateElementDialog,
    CreateSolidElementDialog, CreateShearElementDialog, CreateGapElementDialog,
    CreateBushElementDialog, CreateRbeDialog, CreateConm2Dialog, CreateCoordDialog,
    CreatePlotelDialog,
    CreateSpiderDialog, CreateWeldDialog,
)
from node_runner.dialogs.selection import EntitySelectionDialog
from node_runner.dialogs.tools import (
    CoincidentNodeDialog, GeneratorDialog, ColorManagerDialog,
    NodeTransformDialog, ImportOptionsDialog, DuplicateElementDialog,
    FindReplaceDialog, RenumberDialog, ImportCADDialog,
)
from node_runner.dialogs.info import (
    InfoDialog, LenientImportReportDialog, MaterialImportReportDialog,
    FreeBodyDiagramDialog,
    MassPropertiesReportDialog, FreeEdgeReportDialog,
    QualitySummaryReportDialog, OrphanCheckDialog,
)
from node_runner.dialogs.loads import (
    CreateLoadDialog, CreateConstraintDialog, CreateLoadCombinationDialog,
    SubcaseEditorDialog, CreateEigrlDialog,
)
from node_runner.dialogs.analysis import (
    AnalysisSet, AnalysisSetManagerDialog, OutputRequestsDialog,
    SubcaseMatrixDialog,
)
from node_runner.dialogs.widgets import OrientationWidget
from node_runner.dialogs.geometry import (
    CreateGeometryPointDialog, CreateGeometryLineDialog,
    CreateGeometryArcDialog, CreateGeometryCircleDialog,
    CreateGeometrySurfaceDialog, MeshCurveDialog, MeshSurfaceDialog,
)
from node_runner.dialogs.export import (
    ExportOptionsDialog, ExportDefaultsDialog, UnitsDialog, SaveBdfDialog,
)
from node_runner.dialogs.command_palette import (
    CommandPaletteDialog, install_command_palette_shortcut,
)
from node_runner.dialogs.cross_section import (
    CrossSectionDialog, CrossSectionDock, CrossSectionPanel,
    DefinePlaneDialog,
)
from node_runner.dialogs.mesh_edit import (
    SmoothNodesDialog, MirrorElementsDialog, CopyElementsDialog,
    CombineTriasDialog, InsertEdgeNodeDialog,
)
from node_runner.dialogs.loads_extra import (
    CreatePload1Dialog, CreatePload2Dialog, CreateSpcdDialog,
    CreateMpcDialog, CreateRbarDialog, CreateRbe1Dialog,
    CreateRsplineDialog, CreateBoltDialog,
)
from node_runner.dialogs.result_browser import (
    ResultBrowserDock, AnimationTimelineWidget, VectorOverlayWidget,
)
from node_runner.dialogs.units_conversion import UnitConversionDialog
