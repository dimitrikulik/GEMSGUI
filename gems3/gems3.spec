%define	name	gems
%define	version	3.0.0
%define	release	1497
%define install_dir /usr
%define shared_dir %{install_dir}/share
%define shared_src_dir shared

#%define libname_orig %mklibname %{name}
#%define libname %{libname_orig}0

Name:		%{name}
Summary: 	Chemistry, geochemistry, Gibbs energy minimization
Version:	%{version}
Release:	%{release}
#Copyright:	Distributable
License:	Proprietary_non_commercial
Group:		Applications/Engineering
Source0: 	%{name}-%{version}.tar.bz2
#BuildRoot:	%{_tmppath}/%{name}-%{version}-%{release}-buildroot
#Target: i586-gnu-linux
Url:		http://gems.web.psi.ch/
#BuildRequires: 	libqt3-devel
#Requires:	%{libname} = %{version}-%{release}

%description
This program calculates complex geochemical equilibria
and helps in preparing consistent thermodynamic databases.
The code is developed jointly at PSI LES and ETHZ IGMR,
Switzerland. The code is available free of charge for 
non-profit academic research.
%prep

%setup

%build
#make clean
qmake gems_qt.pro
make
#make RELEASE=YES

%install

#if exists 
/usr/sbin/groupadd gems_admin || /bin/true

install -s -m 755 -o 0 -g 0 gems3 %{install_dir}/bin/gems3

install -d -m 775 -o 0 -g gems_admin %{shared_dir}/gems3
install -d -m 775 -o 0 -g gems_admin %{shared_dir}/gems3/DB.default
install -m 664 -o 0 -g gems_admin %{shared_src_dir}/DB.default/* %{shared_dir}/gems3/DB.default 
install -d -m 775 -o 0 -g gems_admin %{shared_dir}/gems3/projects
install -m 664 -o 0 -g gems_admin %{shared_src_dir}/projects/*.* %{shared_dir}/gems3/projects 
install -d -m 775 -o 0 -g gems_admin %{shared_dir}/gems3/projects/TryNPTDB
install -m 664 -o 0 -g gems_admin %{shared_src_dir}/projects/TryNPTDB/* %{shared_dir}/gems3/projects/TryNPTDB
# Other test/tutorial projects to be added here 
#
install -d -m 775 -o 0 %{shared_dir}/gems3/visor.data
install -m 664 -o 0 %{shared_src_dir}/visor.data/* %{shared_dir}/gems3/visor.data
install -d -m 775 -o 0 %{shared_dir}/gems3/data
install -m 664 -o 0 %{shared_src_dir}/data/* %{shared_dir}/gems3/data
install -d -m 775 -o 0 %{shared_dir}/gems3/doc/html
install -m 664 -o 0 %{shared_src_dir}/doc/html/* %{shared_dir}/gems3/doc/html
install -d -m 775 -o 0 %{shared_dir}/gems3/doc/pdf
install -m 664 -o 0 %{shared_src_dir}/doc/pdf/* %{shared_dir}/gems3/doc/pdf
install -d -m 775 -o 0 %{shared_dir}/gems3/doc/txt
install -m 664 -o 0 %{shared_src_dir}/doc/txt/* %{shared_dir}/gems3/doc/txt

install -d -m 775 -o 0 %{shared_dir}/gems3/img
install -m 664 -o 0 %{shared_src_dir}/img/* %{shared_dir}/gems3/img
install -m 664 -o 0 %{shared_src_dir}/img/gems32.xpm %{shared_dir}/icons/gems3.xpm
install -m 664 -o 0 %{shared_src_dir}/img/gems32.xpm %{shared_dir}/pixmaps/gems3.xpm
#install -m 664 -o 0 img/gems16.xpm %{shared_dir}/icons/mini/gems3.xpm
install -m 664 -o 0 doc/gems3.desktop %{shared_dir}/applications/

%files
%{install_dir}/bin/gems3
%{shared_dir}/gems3
%{shared_dir}/icons/gems3.xpm
%{shared_dir}/pixmaps/gems3.xpm
%{shared_dir}/applications/gems3.desktop