function centerApp(app,parentApp)
    parentPos = parentApp.UIFigure.Position;
    pos = app.UIFigure.Position;
    movegui(app.UIFigure,[round(parentPos(1)+parentPos(3)/2-pos(3)/2),round(parentPos(2)+parentPos(4)/2)-pos(4)/2]);
end