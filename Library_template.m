function MSLib = TBDMS()

MSLib = struct('Name',{},'Derivatization',{},'RetTime',{},...
    'MainIons',{},'MainIonsRelInt',{},'SelectedIons',{},'SelectedIonsFormula',{},'SelectedIonsAtoms',{});

% Fragment name
Item.Name = 'Frag_'; %Enter name of fragment
Item.Derivatization = 'TBDMS';
Item.RetTime = 'retention time'; %enter retention time here
PyrRT = Item.RetTime;
Item.MainIons = ['m0','mf']; %enter m/z range from m0 to mf here
Item.MainIonsRelInt = [100];
Item.SelectedIons = ['m0','mf']; %enter m/z range from m0 to mf here
Item.SelectedIonsFormula = {'formula'}; %enter fragment chemical formula here. For example, pyruvate 174 m/z entered as c6h12o3n1si
Item.SelectedIonsAtoms = {''};
MSLib(end+1) = Item;
