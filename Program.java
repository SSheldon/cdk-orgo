import org.openscience.cdk.Molecule;
import org.openscience.cdk.templates.MoleculeFactory;
//import org.openscience.cdk.io.CMLReader;
//import javax.swing.JFileChooser;
//import java.io.*;

public class Program
{
    public static void main(String[] args)
    {
        OrganicMolecule x = new OrganicMolecule(MoleculeFactory.makeAlkane(7));
        System.out.println(x.getIUPACName());
        //JFileChooser chooser = new JFileChooser();
        //if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
        //{
        //    try
        //    {
        //        CMLReader reader = new CMLReader(new FileInputStream(chooser.getSelectedFile()));
        //    }
        //    catch (FileNotFoundException e) { }
        //}
    }
}
