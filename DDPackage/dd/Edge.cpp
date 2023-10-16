#include "dd/Edge.hpp"

#include "dd/Node.hpp"

namespace dd {

	template <class Node>
	const Edge<Node> Edge<Node>::zero{Node::getTerminal(), Complex::zero, the_maps::the_maps_header()};
	template <class Node>
	const Edge<Node> Edge<Node>::one{Node::getTerminal(), Complex::one, the_maps::the_maps_header()};

	template <class Node> Edge<Node> Edge<Node>::terminal(const Complex& w) {
		return { Node::getTerminal(), w ,the_maps::the_maps_header()};
	}

	template <class Node> bool Edge<Node>::isTerminal() const {
		return Node::isTerminal(p);
	}

	template <class Node> bool Edge<Node>::isZeroTerminal() const {
		return isTerminal() && w == Complex::zero;
	}

	template <class Node> bool Edge<Node>::isOneTerminal() const {
		return isTerminal() && w == Complex::one;
	}

	// Explicit instantiations

	template struct Edge<mNode>;


} // namespace dd


