

function orthoview(img);
clf;
x = round(size(img,1)/2);
y = round(size(img,2)/2);
z = round(size(img,3)/2);

subplot(2,2,2);
hA = imshow(reshape(img(:,:,z),[size(img,1) size(img,2)]),'InitialMagnification',100);
set(hA,'ButtonDownFcn',{@axisA_callback,gcbo,img});
set(hA,'Tag','imageA');
line([1 1 ],[ size(img,1) 0 ],'Color',[0 1 1],'Tag','crosshairAX','ButtonDownFcn',{@line_buttondown,[],img});
line([ size(img,2) 0 ], [1 1],'Color',[0 1 1],'Tag','crosshairAY','ButtonDownFcn',{@line_buttondown,[],img});
title('XY');

subplot(2,2,1);
hB = imshow(reshape(img(:,y,:),[size(img,1) size(img,3)]),'InitialMagnification',100);
set(hB,'ButtonDownFcn',{@axisB_callback,gcbo,img});
set(hB,'Tag','imageB');
line([1 1 ],[ size(img,1) 0 ],'Color',[0 1 1],'Tag','crosshairBX','ButtonDownFcn',{@line_buttondown,[],img});
line([ size(img,2) 0 ], [1 1],'Color',[0 1 1],'Tag','crosshairBY','ButtonDownFcn',{@line_buttondown,[],img});
title('XZ');

subplot(2,2,3);
hC = imshow(reshape(img(x,:,:),[size(img,2) size(img,3)]),'InitialMagnification',100);
set(hC,'ButtonDownFcn',{@axisC_callback,gcbo,img});
set(hC,'Tag','imageC');
line([1 1 ],[ size(img,1) 0 ],'Color',[0 1 1],'Tag','crosshairCX','ButtonDownFcn',{@line_buttondown,[],img});
line([ size(img,2) 0 ], [1 1],'Color',[0 1 1],'Tag','crosshairCY','ButtonDownFcn',{@line_buttondown,[],img});
title('YZ');

subplot(2,2,4);
axis off;
text(0,1,sprintf('dimensions (%i,%i,%i)',size(img,1),size(img,2),size(img,3)));
text(0,0.8,sprintf('%i',x),'Tag','Xinfo');
text(0.2,0.8,sprintf('%i',y),'Tag','Yinfo');
text(0.4,0.8,sprintf('%i',z),'Tag','Zinfo');
text(0,0.6,'x','Tag','Cinfo');


hf = get(gcf, 'javaframe');
rootpane = getAxisComponent(hf);
set(rootpane, 'MouseWheelMovedCallback', {@wheelfcn, gcf} );



updateCinfo(img);
updateAX(y,img);
updateAY(x,img);
updateBX(z,img);
updateBY(x,img);
updateCX(z,img);
updateCY(y,img);



return;


function updateCinfo(img)
hx = str2num(get(findobj(gcf,'Tag','Xinfo'),'String'));
hy = str2num(get(findobj(gcf,'Tag','Yinfo'),'String'));
hz = str2num(get(findobj(gcf,'Tag','Zinfo'),'String'));
set(findobj(gcf,'Tag','Cinfo'),'String',sprintf('Value: %.2f',img(hx,hy,hz)));



function varargout = axisA_callback(h, eventdata, handles, varargin)
pos = get(get(h,'Parent'),'CurrentPoint');
img = varargin{1};
updateAX(pos(1,1),img);
updateAY(pos(1,2),img);
updateCinfo(img);
function varargout = axisB_callback(h, eventdata, handles, varargin)
pos = get(get(h,'Parent'),'CurrentPoint');
img = varargin{1};
updateBX(pos(1,1),img);
updateBY(pos(1,2),img);
updateCinfo(img);
function varargout = axisC_callback(h, eventdata, handles, varargin)
pos = get(get(h,'Parent'),'CurrentPoint');
img = varargin{1};
updateCX(pos(1,1),img);
updateCY(pos(1,2),img);
updateCinfo(img);




function setCrossHair(name,X,Y)
h = findobj(gca,'Tag',name); 
set(h,'XData',X); set(h,'YData',Y);

function varargout = line_buttondown(h, eventdata, handles, varargin)
if length(varargin) > 0,
    set(gcf,'WindowButtonMotionFcn',{@online_mousemove,[],get(h,'Tag'),varargin{1}});
    set(gcf,'WindowButtonUpFcn',{@online_buttonup});
end;

function varargout = online_buttonup(h, eventdata, handles, varargin)
set(gcf,'WindowButtonMotionFcn',{});


function varargout = online_mousemove(h, eventdata, handles, varargin)
if sum(varargin{1} ~= 'crosshairAX') == 0
    pos = get(get(findobj(gcf,'Tag','imageA'),'Parent'),'CurrentPoint');
    updateAX(pos(1,1),varargin{2});
end;
if sum(varargin{1} ~= 'crosshairAY') == 0
    pos = get(get(findobj(gcf,'Tag','imageA'),'Parent'),'CurrentPoint');
    updateAY(pos(1,2),varargin{2});
end;
if sum(varargin{1} ~= 'crosshairBX') == 0
    pos = get(get(findobj(gcf,'Tag','imageB'),'Parent'),'CurrentPoint');
    updateBX(pos(1,1),varargin{2});
end;
if sum(varargin{1} ~= 'crosshairBY') == 0
    pos = get(get(findobj(gcf,'Tag','imageB'),'Parent'),'CurrentPoint');
    updateBY(pos(1,2),varargin{2});
end;
if sum(varargin{1} ~= 'crosshairCX') == 0
    pos = get(get(findobj(gcf,'Tag','imageC'),'Parent'),'CurrentPoint');
    updateCX(pos(1,1),varargin{2});
end;
if sum(varargin{1} ~= 'crosshairCY') == 0
    pos = get(get(findobj(gcf,'Tag','imageC'),'Parent'),'CurrentPoint');
    updateCY(pos(1,2),varargin{2});
end;
updateCinfo(varargin{2});





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sliceupdates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateAX(px,img)
if px < 1 || px > size(img,2),
    return;
end;
subplot(2,2,2);
setCrossHair('crosshairAX',[px px],[ size(img,1) 0 ]);
hy = findobj(gcf,'Tag','Yinfo'); set(hy,'String',sprintf('%i,',round(px)));
subplot(2,2,1);
hB = findobj(gcf,'Tag','imageB');
set(hB,'CData',reshape(img(:,round(px),:),[size(img,1) size(img,3)]));
subplot(2,2,3);
setCrossHair('crosshairCY',[ size(img,3) 0 ],[px px]);

function updateAY(py,img)
if py < 1 || py > size(img,1),
    return;
end;
subplot(2,2,2);
setCrossHair('crosshairAY',[ size(img,2) 0 ],[py py]);
hx = findobj(gcf,'Tag','Xinfo'); set(hx,'String',sprintf('%i,',round(py)));
subplot(2,2,1);
setCrossHair('crosshairBY',[ size(img,3) 0 ],[py py]);
subplot(2,2,3);
hC = findobj(gcf,'Tag','imageC');
set(hC,'CData',reshape(img(round(py),:,:),[size(img,2) size(img,3)]));

function updateBX(px,img)
if px < 1 || px > size(img,3),
    return;
end;
subplot(2,2,1);
hz = findobj(gcf,'Tag','Zinfo'); set(hz,'String',sprintf('%i',round(px)));
setCrossHair('crosshairBX',[px px ],[ size(img,1) 0 ]);
subplot(2,2,2);
hA = findobj(gcf,'Tag','imageA');
set(hA,'CData',reshape(img(:,:,round(px)),[size(img,1) size(img,2)]));
subplot(2,2,3);
setCrossHair('crosshairCX',[px px],[ size(img,3) 0 ]);

function updateBY(py,img)
if py < 1 || py > size(img,1),
    return;
end;
subplot(2,2,1);
hx = findobj(gcf,'Tag','Xinfo'); set(hx,'String',sprintf('%i,',round(py)));
setCrossHair('crosshairBY',[ size(img,3) 0 ], [py py]);
subplot(2,2,2);
setCrossHair('crosshairAY',[ size(img,1) 0 ],[py py ]);
subplot(2,2,3);
hC = findobj(gcf,'Tag','imageC');
set(hC,'CData',reshape(img(round(py),:,:),[size(img,2) size(img,3)]));

function updateCX(px,img)
if px < 1 || px > size(img,3),
    return;
end;
subplot(2,2,3);
setCrossHair('crosshairCX',[px px],[ size(img,1) 0 ]);
hz = findobj(gcf,'Tag','Zinfo'); set(hz,'String',sprintf('%i',round(px)));
subplot(2,2,2);
hA = findobj(gcf,'Tag','imageA');
set(hA,'CData',reshape(img(:,:,round(px)),[size(img,1) size(img,2)]));
subplot(2,2,1);
setCrossHair('crosshairBX',[px px],[ size(img,3) 0 ]);

function updateCY(py,img)
if py < 1 || py > size(img,2),
    return;
end;
subplot(2,2,3);
setCrossHair('crosshairCY',[ size(img,3) 0 ], [py py]);
hy = findobj(gcf,'Tag','Yinfo'); set(hy,'String',sprintf('%i,',round(py)));
subplot(2,2,2);
setCrossHair('crosshairAX',[py py],[ size(img,1) 0 ]);
subplot(2,2,1);
hB = findobj(gcf,'Tag','imageB');
set(hB,'CData',reshape(img(:,round(py),:),[size(img,1) size(img,3)]));









function wheelfcn( dummy, eventdata, hObj ) %#ok
eventdata = struct(get(eventdata));
pw_pos = getpixelposition(gcf); %get pixel position of this figure
pl = get(0,'PointerLocation'); %get screen location of pointer
point = pl - pw_pos(1:2); %get pointer location relative to window
imH = hittest(gcf,point);
hx = str2num(get(findobj(gcf,'Tag','Xinfo'),'String'));
hy = str2num(get(findobj(gcf,'Tag','Yinfo'),'String'));
hz = str2num(get(findobj(gcf,'Tag','Zinfo'),'String'));
dummy = get(findobj(gcf,'Tag','imageA'),'ButtonDownFcn');
img = dummy{3};

step = eventdata.WheelRotation*5;
if eventdata.Modifiers == 1,
    step = step / 5;
end;
switch get(imH,'Tag'),
    case {'imageA','crosshairAX','crosshairAY'}
        hz = hz + step;
        updateBX(hz,img);
    case {'imageB','crosshairBX','crosshairBY'}
        hy = hy + step;
        updateCY(hy,img);
    case {'imageC','crosshairCX','crosshairCY'}
        hx = hx + step;
        updateAY(hx,img);   
end;








