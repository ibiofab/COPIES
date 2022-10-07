import './App.css';
import OrganismInfo from './components/Organism Info/OrganismInfo';
import HomologyArmInfo from './components/Homology Arm Info/HomologyArmInfo';
import HarbourInfo from './components/Harbour Info/HarbourInfo';
import GuideRNAInfo from './components/Guide RNA Info/GuideRNAInfo';

function App() {
  return (
    <div className="App">
      <OrganismInfo />
      <HomologyArmInfo /> 
      <GuideRNAInfo />
      <HarbourInfo />
    </div>
  );
}

export default App;
