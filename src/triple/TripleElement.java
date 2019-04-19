package triple;

/**
 * An object of this class represents a element that could be inserted in a prohibited list
 * for MAXQBFPT problem. An element is simply a variable of the instance.
 * 
 * @author Cintia Muranaka
 * @author Felipe de Carvalho Pereira [felipe.pereira@students.ic.unicamp.br]
 * @author Matheus Diógenes Andrade
 */
public class TripleElement {

    public final Integer index; //The index of the variable
    public Boolean selected; //If the element is already selected in the partial solution
    public Boolean available; // If the element is available to be inserted into the partial solution
    public Integer incumbentFrequency; // The number of times that this element was present in the incumbent solutions



	public TripleElement(int index) {
        this.index = index;
        this.selected = false;
        this.available = true;
        this.incumbentFrequency = 0;
    }
	
	// Increase frequency in 1
	public void increaseFrequency()
	{
		this.incumbentFrequency++;
	}
	
	//Getters and Setters
	
	public Integer getIndex() {
		return index;
	}
	
	public Boolean getSelected() {
		return selected;
	}

	public Boolean getAvailable() {
		return available;
	}
	
	public Integer getIncumbentFrequency() {
		return incumbentFrequency;
	}

	public void setAvailable(Boolean available) {
		this.available = available;
	}
	
	public void setSelected(Boolean selected) {
		this.selected = selected;
	}
	
    public void setIncumbentFrequency(Integer incumbentFrequency) {
		this.incumbentFrequency = incumbentFrequency;
	}

}
