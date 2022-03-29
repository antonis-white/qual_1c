#include "subway.h"

using namespace subway;
using namespace std;

int main() {
	line_t line1({0, 1, 2}, {2, 2}, 2);
	line_t line2({3, 1, 4}, {3, 3}, 3);
	subway_predictor pred(5, {line1, line2});
	cout << pred.predict(0, 2, 0) << '\n';
	cout << pred.predict(0, 2, 1) << '\n';
	cout << pred.predict(0, 3, 0) << '\n';
	cout << pred.predict(0, 3, 1) << '\n';
}