# -*- coding: utf-8 -*-

# This file is part of pyqt4topyqt5

MODULES = [
    'QAxContainer',
    'QtCore',
    'QtDBus',
    'QtDesigner',
    'QtGui',
    'QtHelp',
    'QtMultimedia',
    'QtNetwork',
    'QtOpenGL',
    'QtPrintSupport',
    'QtSql',
    'QtSvg',
    'QtTest',
    'QtWebKit',
    'QtWebKitWidgets',
    'QtWidgets',
    'QtXmlPatterns',
    ]

DISCARDED = {
    'QtOpenGl': [
        'QGLBuffer',
        'QGLColormap',
        'QGLFramebufferObject',
        'QGLFramebufferObjectFormat',
        'QGLPixelBuffer',
        'QGLShader',
        'QGLShaderProgram']}

CLASSES = {
    'QAxContainer': [
        'QAxBase',
        'QAxObject',
        'QAxWidget'],
    'QtCore': [
        'QAbstractAnimation',
        'QAbstractEventDispatcher',
        'QAbstractItemModel',
        'QAbstractListModel',
        'QAbstractNativeEventFilter',
        'QAbstractProxyModel',
        'QAbstractState',
        'QAbstractTableModel',
        'QAbstractTransition',
        'QAnimationGroup',
        'QBasicTimer',
        'QBitArray',
        'QBuffer',
        'QByteArray',
        'QByteArrayMatcher',
        'QChildEvent',
        'QCoreApplication',
        'QCryptographicHash',
        'QDataStream',
        'QDate',
        'QDateTime',
        'QDir',
        'QDirIterator',
        'QDynamicPropertyChangeEvent',
        'QEasingCurve',
        'QElapsedTimer',
        'QEvent',
        'QEventLoop',
        'QEventLoopLocker',
        'QEventTransition',
        'QFile',
        'QFileDevice',
        'QFileInfo',
        'QFileSystemWatcher',
        'QFinalState',
        'QGenericArgument',
        'QGenericReturnArgument',
        'QHistoryState',
        'QIODevice',
        'QIdentityProxyModel',
        'QItemSelection',
        'QItemSelectionModel',
        'QItemSelectionRange',
        'QLibrary',
        'QLibraryInfo',
        'QLine',
        'QLineF',
        'QLocale',
        'QMargins',
        'QMessageLogContext',
        'QMessageLogger',
        'QMetaClassInfo',
        'QMetaEnum',
        'QMetaMethod',
        'QMetaObject',
        'QMetaProperty',
        'QMetaType',
        'QMimeData',
        'QMimeDatabase',
        'QMimeType',
        'QModelIndex',
        'QMutex',
        'QMutexLocker',
        'QObject',
        'QObjectCleanupHandler',
        'QParallelAnimationGroup',
        'QPauseAnimation',
        'QPersistentModelIndex',
        'QPluginLoader',
        'QPoint',
        'QPointF',
        'QProcess',
        'QProcessEnvironment',
        'QPropertyAnimation',
        'QReadLocker',
        'QReadWriteLock',
        'QRect',
        'QRectF',
        'QRegExp',
        'QRegularExpression',
        'QRegularExpressionMatch',
        'QRegularExpressionMatchIterator',
        'QResource',
        'QRunnable',
        'QSemaphore',
        'QSequentialAnimationGroup',
        'QSettings',
        'QSharedMemory',
        'QSignalMapper',
        'QSignalTransition',
        'QSize',
        'QSizeF',
        'QSocketNotifier',
        'QSortFilterProxyModel',
        'QStandardPaths',
        'QState',
        'QStateMachine',
        'QStringListModel',
        'QSysInfo',
        'QSystemSemaphore',
        'QTemporaryDir',
        'QTemporaryFile',
        'QTextBoundaryFinder',
        'QTextCodec',
        'QTextDecoder',
        'QTextEncoder',
        'QTextStream',
        'QTextStreamManipulator',
        'QThread',
        'QThreadPool',
        'QTime',
        'QTimeLine',
        'QTimer',
        'QTimerEvent',
        'QTranslator',
        'QUrl',
        'QUrlQuery',
        'QUuid',
        'QVariant',
        'QVariantAnimation',
        'QWaitCondition',
        'QWriteLocker',
        'QXmlStreamAttribute',
        'QXmlStreamAttributes',
        'QXmlStreamEntityDeclaration',
        'QXmlStreamEntityResolver',
        'QXmlStreamNamespaceDeclaration',
        'QXmlStreamNotationDeclaration',
        'QXmlStreamReader',
        'QXmlStreamWriter',
        'Qt'],
    'QtDBus': [
        'QDBus',
        'QDBusAbstractAdaptor',
        'QDBusAbstractInterface',
        'QDBusArgument',
        'QDBusConnection',
        'QDBusConnectionInterface',
        'QDBusError',
        'QDBusInterface',
        'QDBusMessage',
        'QDBusObjectPath',
        'QDBusPendingCall',
        'QDBusPendingCallWatcher',
        'QDBusPendingReply',
        'QDBusReply',
        'QDBusServiceWatcher',
        'QDBusSignature',
        'QDBusUnixFileDescriptor',
        'QDBusVariant'],
    'QtDesigner': [
        'QAbstractExtensionFactory',
        'QAbstractExtensionManager',
        'QAbstractFormBuilder',
        'QDesignerActionEditorInterface',
        'QDesignerContainerExtension',
        'QDesignerCustomWidgetCollectionInterface',
        'QDesignerCustomWidgetInterface',
        'QDesignerFormEditorInterface',
        'QDesignerFormWindowCursorInterface',
        'QDesignerFormWindowInterface',
        'QDesignerFormWindowManagerInterface',
        'QDesignerMemberSheetExtension',
        'QDesignerObjectInspectorInterface',
        'QDesignerPropertyEditorInterface',
        'QDesignerPropertySheetExtension',
        'QDesignerTaskMenuExtension',
        'QDesignerWidgetBoxInterface',
        'QExtensionFactory',
        'QExtensionManager',
        'QFormBuilder'],
    'QtGui': [
        'QAbstractTextDocumentLayout',
        'QActionEvent',
        'QBackingStore',
        'QBitmap',
        'QBrush',
        'QClipboard',
        'QCloseEvent',
        'QColor',
        'QConicalGradient',
        'QContextMenuEvent',
        'QCursor',
        'QDesktopServices',
        'QDoubleValidator',
        'QDrag',
        'QDragEnterEvent',
        'QDragLeaveEvent',
        'QDragMoveEvent',
        'QDropEvent',
        'QEnterEvent',
        'QExposeEvent',
        'QFileOpenEvent',
        'QFocusEvent',
        'QFont',
        'QFontDatabase',
        'QFontInfo',
        'QFontMetrics',
        'QFontMetricsF',
        'QGlyphRun',
        'QGradient',
        'QGuiApplication',
        'QHelpEvent',
        'QHideEvent',
        'QHoverEvent',
        'QIcon',
        'QIconDragEvent',
        'QIconEngine',
        'QImage',
        'QImageIOHandler',
        'QImageReader',
        'QImageWriter',
        'QInputEvent',
        'QInputMethod',
        'QInputMethodEvent',
        'QInputMethodQueryEvent',
        'QIntValidator',
        'QKeyEvent',
        'QKeySequence',
        'QLinearGradient',
        'QMatrix2x2',
        'QMatrix2x3',
        'QMatrix2x4',
        'QMatrix3x2',
        'QMatrix3x3',
        'QMatrix3x4',
        'QMatrix4x2',
        'QMatrix4x3',
        'QMatrix4x4',
        'QMouseEvent',
        'QMoveEvent',
        'QMovie',
        'QOpenGLContext',
        'QOpenGLContextGroup',
        'QPagedPaintDevice',
        'QPaintDevice',
        'QPaintEngine',
        'QPaintEngineState',
        'QPaintEvent',
        'QPainter',
        'QPainterPath',
        'QPainterPathStroker',
        'QPalette',
        'QPdfWriter',
        'QPen',
        'QPicture',
        'QPictureIO',
        'QPixmap',
        'QPixmapCache',
        'QPolygon',
        'QPolygonF',
        'QPyTextObject',
        'QQuaternion',
        'QRadialGradient',
        'QRawFont',
        'QRegExpValidator',
        'QRegion',
        'QResizeEvent',
        'QScreen',
        'QScrollEvent',
        'QScrollPrepareEvent',
        'QSessionManager',
        'QShortcutEvent',
        'QShowEvent',
        'QStandardItem',
        'QStandardItemModel',
        'QStaticText',
        'QStatusTipEvent',
        'QStyleHints',
        'QSurface',
        'QSurfaceFormat',
        'QSyntaxHighlighter',
        'QTabletEvent',
        'QTextBlock',
        'QTextBlockFormat',
        'QTextBlockGroup',
        'QTextBlockUserData',
        'QTextCharFormat',
        'QTextCursor',
        'QTextDocument',
        'QTextDocumentFragment',
        'QTextDocumentWriter',
        'QTextFormat',
        'QTextFragment',
        'QTextFrame',
        'QTextFrameFormat',
        'QTextImageFormat',
        'QTextInlineObject',
        'QTextItem',
        'QTextLayout',
        'QTextLength',
        'QTextLine',
        'QTextList',
        'QTextListFormat',
        'QTextObject',
        'QTextObjectInterface',
        'QTextOption',
        'QTextTable',
        'QTextTableCell',
        'QTextTableCellFormat',
        'QTextTableFormat',
        'QTouchDevice',
        'QTouchEvent',
        'QTransform',
        'QValidator',
        'QVector2D',
        'QVector3D',
        'QVector4D',
        'QWhatsThisClickedEvent',
        'QWheelEvent',
        'QWindow',
        'QWindowStateChangeEvent'],
    'QtHelp': [
        'QHelpContentItem',
        'QHelpContentModel',
        'QHelpContentWidget',
        'QHelpEngine',
        'QHelpEngineCore',
        'QHelpIndexModel',
        'QHelpIndexWidget',
        'QHelpSearchEngine',
        'QHelpSearchQuery',
        'QHelpSearchQueryWidget',
        'QHelpSearchResultWidget'],
    'QtMultimedia': [
        'QAbstractVideoBuffer',
        'QAbstractVideoSurface',
        'QAudio',
        'QAudioDeviceInfo',
        'QAudioFormat',
        'QAudioInput',
        'QAudioOutput',
        'QVideoFrame',
        'QVideoSurfaceFormat'],
    'QtNetwork': [
        'QAbstractNetworkCache',
        'QAbstractSocket',
        'QAuthenticator',
        'QDnsDomainNameRecord',
        'QDnsHostAddressRecord',
        'QDnsLookup',
        'QDnsMailExchangeRecord',
        'QDnsServiceRecord',
        'QDnsTextRecord',
        'QHostAddress',
        'QHostInfo',
        'QHttpMultiPart',
        'QHttpPart',
        'QLocalServer',
        'QLocalSocket',
        'QNetworkAccessManager',
        'QNetworkAddressEntry',
        'QNetworkCacheMetaData',
        'QNetworkConfiguration',
        'QNetworkConfigurationManager',
        'QNetworkCookie',
        'QNetworkCookieJar',
        'QNetworkDiskCache',
        'QNetworkInterface',
        'QNetworkProxy',
        'QNetworkProxyFactory',
        'QNetworkProxyQuery',
        'QNetworkReply',
        'QNetworkRequest',
        'QNetworkSession',
        'QSsl',
        'QSslCertificate',
        'QSslCertificateExtension',
        'QSslCipher',
        'QSslConfiguration',
        'QSslError',
        'QSslKey',
        'QSslSocket',
        'QTcpServer',
        'QTcpSocket',
        'QUdpSocket'],
    'QtOpenGL': [
        'QGL',
        'QGLContext',
        'QGLFormat',
        'QGLWidget'],
    'QtPrintSupport': [
        'QAbstractPrintDialog',
        'QPageSetupDialog',
        'QPrintDialog',
        'QPrintEngine',
        'QPrinter',
        'QPrinterInfo',
        'QPrintPreviewDialog',
        'QPrintPreviewWidget'],
    'QtSql': [
        'QSql',
        'QSqlDatabase',
        'QSqlDriver',
        'QSqlDriverCreatorBase',
        'QSqlError',
        'QSqlField',
        'QSqlIndex',
        'QSqlQuery',
        'QSqlQueryModel',
        'QSqlRecord',
        'QSqlRelation',
        'QSqlRelationalDelegate',
        'QSqlRelationalTableModel',
        'QSqlResult',
        'QSqlTableModel'],
    'QtSvg': [
        'QGraphicsSvgItem',
        'QSvgGenerator',
        'QSvgRenderer',
        'QSvgWidget'],
    'QtTest': [
        'QSignalSpy',
        'QTest'],
    'QtWebKit': [
        'QWebDatabase',
        'QWebElement',
        'QWebElementCollection',
        'QWebHistory',
        'QWebHistoryInterface',
        'QWebHistoryItem',
        'QWebPluginFactory',
        'QWebSecurityOrigin',
        'QWebSettings'],
    'QtWebKitWidgets': [
        'QGraphicsWebView',
        'QWebFrame',
        'QWebHitTestResult',
        'QWebInspector',
        'QWebPage',
        'QWebView'],
    'QtWidgets': [
        'QAbstractButton',
        'QAbstractGraphicsShapeItem',
        'QAbstractItemDelegate',
        'QAbstractItemView',
        'QAbstractScrollArea',
        'QAbstractSlider',
        'QAbstractSpinBox',
        'QAction',
        'QActionGroup',
        'QApplication',
        'qApp',
        'QApplication.instance()',
        'QBoxLayout',
        'QButtonGroup',
        'QCalendarWidget',
        'QCheckBox',
        'QColorDialog',
        'QColumnView',
        'QComboBox',
        'QCommandLinkButton',
        'QCommonStyle',
        'QCompleter',
        'QDataWidgetMapper',
        'QDateEdit',
        'QDateTimeEdit',
        'QDesktopWidget',
        'QDial',
        'QDialog',
        'QDialogButtonBox',
        'QDirModel',
        'QDockWidget',
        'QDoubleSpinBox',
        'QErrorMessage',
        'QFileDialog',
        'QFileIconProvider',
        'QFileSystemModel',
        'QFocusFrame',
        'QFontComboBox',
        'QFontDialog',
        'QFormLayout',
        'QFrame',
        'QGesture',
        'QGestureEvent',
        'QGestureRecognizer',
        'QGraphicsAnchor',
        'QGraphicsAnchorLayout',
        'QGraphicsBlurEffect',
        'QGraphicsColorizeEffect',
        'QGraphicsDropShadowEffect',
        'QGraphicsEffect',
        'QGraphicsEllipseItem',
        'QGraphicsGridLayout',
        'QGraphicsItem',
        'QGraphicsItemGroup',
        'QGraphicsLayout',
        'QGraphicsLayoutItem',
        'QGraphicsLineItem',
        'QGraphicsLinearLayout',
        'QGraphicsObject',
        'QGraphicsOpacityEffect',
        'QGraphicsPathItem',
        'QGraphicsPixmapItem',
        'QGraphicsPolygonItem',
        'QGraphicsProxyWidget',
        'QGraphicsRectItem',
        'QGraphicsRotation',
        'QGraphicsScale',
        'QGraphicsScene',
        'QGraphicsSceneContextMenuEvent',
        'QGraphicsSceneDragDropEvent',
        'QGraphicsSceneEvent',
        'QGraphicsSceneHelpEvent',
        'QGraphicsSceneHoverEvent',
        'QGraphicsSceneMouseEvent',
        'QGraphicsSceneMoveEvent',
        'QGraphicsSceneResizeEvent',
        'QGraphicsSceneWheelEvent',
        'QGraphicsSimpleTextItem',
        'QGraphicsTextItem',
        'QGraphicsTransform',
        'QGraphicsView',
        'QGraphicsWidget',
        'QGridLayout',
        'QGroupBox',
        'QHBoxLayout',
        'QHeaderView',
        'QInputDialog',
        'QItemDelegate',
        'QItemEditorCreatorBase',
        'QItemEditorFactory',
        'QKeyEventTransition',
        'QLCDNumber',
        'QLabel',
        'QLayout',
        'QLayoutItem',
        'QLineEdit',
        'QListView',
        'QListWidget',
        'QListWidgetItem',
        'QMainWindow',
        'QMdiArea',
        'QMdiSubWindow',
        'QMenu',
        'QMenuBar',
        'QMessageBox',
        'QMouseEventTransition',
        'QPanGesture',
        'QPinchGesture',
        'QPlainTextDocumentLayout',
        'QPlainTextEdit',
        'QProgressBar',
        'QProgressDialog',
        'QPushButton',
        'QRadioButton',
        'QRubberBand',
        'QScrollArea',
        'QScrollBar',
        'QScroller',
        'QScrollerProperties',
        'QShortcut',
        'QSizeGrip',
        'QSizePolicy',
        'QSlider',
        'QSpacerItem',
        'QSpinBox',
        'QSplashScreen',
        'QSplitter',
        'QSplitterHandle',
        'QStackedLayout',
        'QStackedWidget',
        'QStatusBar',
        'QStyle',
        'QStyleFactory',
        'QStyleHintReturn',
        'QStyleHintReturnMask',
        'QStyleHintReturnVariant',
        'QStyleOption',
        'QStyleOptionButton',
        'QStyleOptionComboBox',
        'QStyleOptionComplex',
        'QStyleOptionDockWidget',
        'QStyleOptionFocusRect',
        'QStyleOptionFrame',
        'QStyleOptionGraphicsItem',
        'QStyleOptionGroupBox',
        'QStyleOptionHeader',
        'QStyleOptionMenuItem',
        'QStyleOptionProgressBar',
        'QStyleOptionRubberBand',
        'QStyleOptionSizeGrip',
        'QStyleOptionSlider',
        'QStyleOptionSpinBox',
        'QStyleOptionTab',
        'QStyleOptionTabBarBase',
        'QStyleOptionTabWidgetFrame',
        'QStyleOptionTitleBar',
        'QStyleOptionToolBar',
        'QStyleOptionToolBox',
        'QStyleOptionToolButton',
        'QStyleOptionViewItem',
        'QStylePainter',
        'QStyledItemDelegate',
        'QSwipeGesture',
        'QSystemTrayIcon',
        'QTabBar',
        'QTabWidget',
        'QTableView',
        'QTableWidget',
        'QTableWidgetItem',
        'QTableWidgetSelectionRange',
        'QTapAndHoldGesture',
        'QTapGesture',
        'QTextBrowser',
        'QTextEdit',
        'QTimeEdit',
        'QToolBar',
        'QToolBox',
        'QToolButton',
        'QToolTip',
        'QTreeView',
        'QTreeWidget',
        'QTreeWidgetItem',
        'QTreeWidgetItemIterator',
        'QUndoCommand',
        'QUndoGroup',
        'QUndoStack',
        'QUndoView',
        'QVBoxLayout',
        'QWhatsThis',
        'QWidget',
        'QWidgetAction',
        'QWidgetItem',
        'QWizard',
        'QWizardPage'],
    'QtXmlPatterns': [
        'QAbstractMessageHandler',
        'QAbstractUriResolver',
        'QAbstractXmlNodeModel',
        'QAbstractXmlReceiver',
        'QSimpleXmlNodeModel',
        'QSourceLocation',
        'QXmlFormatter',
        'QXmlItem',
        'QXmlName',
        'QXmlNamePool',
        'QXmlNodeModelIndex',
        'QXmlQuery',
        'QXmlResultItems',
        'QXmlSchema',
        'QXmlSchemaValidator',
        'QXmlSerializer']}

QAPP_STATIC_METHODS = [
    # QCoreApplication
    'addLibraryPath',
    'applicationDirPath',
    'applicationFilePath',
    'applicationName',
    'applicationPid',
    'applicationVersion',
    'argc',
    'arguments',
    'argv',
    'closingDown',
    'exec_',
    'exit',
    'flush',
    'hasPendingEvents',
    'installTranslator',
    'instance',
    'libraryPaths',
    'organizationDomain',
    'organizationName',
    'postEvent',
    'processEvents',
    'quit',
    'removeLibraryPath',
    'removePostedEvents',
    'removeTranslator',
    'sendEvent',
    'sendPostedEvents',
    'setApplicationName',
    'setApplicationVersion',
    'setAttribute',
    'setLibraryPaths',
    'setOrganizationDomain',
    'setOrganizationName',
    'startingUp',
    'testAttribute',
    'translate',

    # QApplication
    'aboutQt',
    'activeModalWidget',
    'activePopupWidget',
    'activeWindow',
    'alert',
    'allWidgets',
    'beep',
    'changeOverrideCursor',
    'clipboard',
    'closeAllWindows',
    'colorSpec',
    'cursorFlashTime',
    'desktop',
    'desktopSettingsAware',
    'doubleClickInterval',
    'focusWidget',
    'font',
    'fontMetrics',
    'globalStrut',
    'isEffectEnabled',
    'isLeftToRight',
    'isRightToLeft',
    'keyboardInputDirection',
    'keyboardInputInterval',
    'keyboardInputLocale',
    'keyboardModifiers',
    'layoutDirection',
    'mouseButtons',
    'overrideCursor',
    'palette',
    'queryKeyboardModifiers',
    'quitOnLastWindowClosed',
    'restoreOverrideCursor',
    'setActiveWindow',
    'setColorSpec',
    'setCursorFlashTime',
    'setDesktopSettingsAware',
    'setDoubleClickInterval',
    'setEffectEnabled',
    'setFont',
    'setGlobalStrut',
    'setGraphicsSystem',
    'setKeyboardInputInterval',
    'setLayoutDirection',
    'setOverrideCursor',
    'setPalette',
    'setQuitOnLastWindowClosed',
    'setStartDragDistance',
    'setStartDragTime',
    'setStyle',
    'setWheelScrollLines',
    'setWindowIcon',
    'startDragDistance',
    'startDragTime',
    'style',
    'syncX',
    'topLevelAt',
    'topLevelWidgets',
    'type',
    'wheelScrollLines',
    'widgetAt',
    'windowIcon',
    ]

QVARIANT_OBSOLETE_METHODS = [
    'toBitArray',
    'toChar',
    'toEasingCurve',
    'toHash',
    'toLineF',
    'toLocale',
    'toMap',
    'toModelIndex',
    'toPyObject',
    'toReal',
    'toRectF',
    'toRegExp',
    'toSizeF',
    'toUrl',
    'toUuid',
    # Many conversion methods that were used by QVariant are still used by other classes
    # and cannot be removed indiscriminately.
    # 'toBool',         # used by QJasonValue and QJSValue
    # 'toByteArray',    # used by QUuid, QNdefMessage, QDomDocument
    # 'toDate',         # used by QLocale
    # 'toDateTime',     # used by QLocale and QJSValue
    # 'toDouble',       # used by QByteArray, QJasonValue and QLocale
    # 'toFloat',        # used by QByteArray and QLocale
    # 'toInt',          # used by QByteArray, QJasonValue, QLocale and QJSValue
    # 'toLine',         # used by QLine
    # 'toList',         # used by QWebEngineScriptCollection and QWebElementCollection
    # 'toLongLong',     # used by QByteArray and QLocale
    # 'toULongLong',    # used by QByteArray and QLocale
    # 'toPoint',        # used by QPoint and QVector[234]D
    # 'toPointF',       # used by QVector[234]D
    # 'toRect',         # used by QRect
    # 'toSize',         # used by QSize
    # 'toString',       # used by too many to list
    # 'toStringList',   # used by QProcessEnvironment and QUrl
    # 'toTime',         # used by QLocale
    # 'toUInt',         # used by QByteArray, QLocale and QJSValue
    ]
