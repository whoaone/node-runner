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
    ExportOptionsDialog, ExportDefaultsDialog, UnitsDialog,
)
from node_runner.dialogs.command_palette import (
    CommandPaletteDialog, install_command_palette_shortcut,
)
from node_runner.dialogs.cross_section import CrossSectionDialog
from node_runner.dialogs.mesh_edit import (
    SmoothNodesDialog, MirrorElementsDialog, CopyElementsDialog,
    CombineTriasDialog,
)
