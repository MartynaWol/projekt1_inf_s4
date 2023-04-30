Ten program został stworzony do przeprowadzania transformacji współrzędnych pomiędzy różnymi układami. Możliwe jest przeprowadzenie transformacji:

	- XYZ (geocentryczne) -> BLH (elipsoidalne fi, lambda, h)
	- BLH -> XYZ
	- XYZ -> NEUp (topocentryczne northing, easting, up)
	- BL(GRS80, WGS84, Krasowski) -> 2000
	- BL(GRS80, WGS84, Krasowski) -> 1992 

Program został napisany dla systemu operacyjnego Windows 10, w programie Spyder w języku programowania Python. Wymagania, które nalezy spełnić aby program działał prawidłowo:

	- posiadać zainstalowany program, w którym zaimplemenowna jest obsługa języka Python (preferowanie Spyder)
	- posiadać zainstalowany program git
	- posiadać konto na stronie github

Instrukcja użycia programu:
	1. Pobierz kod źródłowy programu z repozyturium z GitHub
	2. Otwórz program git i pobierz folder z gałezi main
	3. W tym samym folderze, w którym masz zapisany program utwórz plik z danymi do transformacji
	4. Otwórz Wiersz poleceń
	5. W wierszu poleceń otwórz odpowiednią ścierzkę (tam gdzie masz zapisany program)
	6. Do wiersza poleceń wpisz:
		Python proj1.py <nazwa pliku z danymi>
	7. Program wykona transformację i zapisze plik z wynikami do folderu, w którym pracujesz

Instrukcja tworzenia pliku z danymi:
	1. Utwórz plik tekstowy
	2. Nazwij go odpowiednio. Aby wykonać transformację:
		- XYZ -> BLH: "dane-XYZ_BLH.txt"
		- BLH -> XYZ: "dane-BLH_XYZ.txt"
		- XYZ -> NEU: "dane-XYZ_NEU.txt"
		- BL -> 2000: "dane-BL_U2000.txt"
		- BL -> 1992: "dane-BL_U1992.txt"
	3. W pliku tekstowym zapisz współrzędne które chcesz przetransformować. W przypadku transformacji z BL do U2000 i do U1992 wprowadź też parametry elipsoidy, na któej chcesz wykonać transformacje.
	Dane mają być ustawione w jednej linijce i oddzielone tabulatorami. Po koleji wprowadź:
		- XYZ -> BLH: X, Y, Z
		- BLH -> XYZ: B, L, H
		- XYZ -> NEU: Xa, Ya, Za, Xb, Yb, Zb
		- BL -> 2000: B, L, a, e2
		- BL -> 1992: B, L, a, e2
	*Wartości kątowe zapisuj w radianach, a resztę w metrach

Przykładowy tekst pliku z danymi dla transformacji XYZ -> BLH:

3782480.000	1084940.000	5003080.000

Przykładowy tekst pliku z wynikami transformacji XYZ -> BLH:

X = 3782480.00001056, Y = 1084940.000003029, Z = 5003080.000013967


Znane błędy i nietypowe zachowania programu:
	- program może liczyć błędnie dla określonych niepoprawnych danych wejściowych podanych przez użytkownika, np. gdy wartości długości lub szerokości geograficznej przekroczą 180 stopni, lub gdy wysokość 
	zostanie podana jako 0
	- w przypadku nieprawidłowego zapisu pliku z danymi program może nie wykonać transformacji, lub zostanie ona wykanana błędnie
	- w przypadku umieszczenia pliku z programem i pliku z danymi w innym folderze program nie wykona transformacji
	- program nie wykonuje transformacji wielu punktów na raz