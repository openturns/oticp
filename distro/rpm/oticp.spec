# norootforbuild
%{!?python_sitelib: %global python_sitelib %(%{__python} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")}

Name:           oticp
Version:        0.1
Release:        0%{?dist}
Summary:        Template OpenTURNS module
Group:          System Environment/Libraries
License:        LGPLv3+
URL:            http://www.openturns.org/
Source0:        http://downloads.sourceforge.net/openturns-modules/oticp/oticp-%{version}.tar.bz2
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires:  python-openturns-devel
BuildRequires:  cmake
Requires:       python-oticp

%description
Extra template module for OpenTURNS

%package -n python-%{name}
Summary:        Foo python module
Group:          Productivity/Scientific/Math
Requires:       python-openturns

%description -n python-%{name}
Python textual interface to OTICP library

%prep
%setup -q

%build
cmake -DCMAKE_INSTALL_PREFIX=%{_prefix} \
      -DBUILD_DOC=FALSE .

%install
rm -rf %{buildroot}
make install DESTDIR=%{buildroot}

%clean
rm -rf %{buildroot}

%files -n python-%{name}
%defattr(-,root,root,-)
%{python_sitelib}/
%{_datadir}/%{name}/

%changelog
* Wed Nov 28 2012 Julien Schueller <schueller at phimeca dot com> 0.0-1
- Initial package creation

