// This file is released under the 3-clause BSD license. See COPYING-BSD.
// Generated by builder.sce : Please, do not edit this file
// ----------------------------------------------------------------------------
//
if win64() then
  warning(_("This module requires a Windows x86 platform."));
  return
end
//
fooc_path = get_absolute_file_path('loader.sce');
//
// ulink previous function with same name
[bOK, ilib] = c_link('fooc');
if bOK then
  ulink(ilib);
end
//
link(fooc_path + 'libfooc' + getdynlibext(), ['fooc'],'c');
// remove temp. variables on stack
clear fooc_path;
clear bOK;
clear ilib;
// ----------------------------------------------------------------------------
