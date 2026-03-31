import { Routes, Route } from "react-router-dom";
import Navbar from "./components/Navbar.jsx";
import Footer from "./components/Footer.jsx";
import Landing from "./pages/Landing.jsx";
import Installation from "./pages/Installation.jsx";
import SetupAssistant from "./pages/SetupAssistant.jsx";
import ConfigBuilder from "./pages/ConfigBuilder.jsx";
import Examples from "./pages/Examples.jsx";
import Simulation from "./pages/Simulation.jsx";

export default function App() {
  return (
    <div className="flex min-h-screen flex-col">
      <Navbar />
      <main className="flex-1">
        <Routes>
          <Route path="/" element={<Landing />} />
          <Route path="/install" element={<Installation />} />
          <Route path="/setup" element={<SetupAssistant />} />
          <Route path="/config" element={<ConfigBuilder />} />
          <Route path="/examples" element={<Examples />} />
          <Route path="/simulation" element={<Simulation />} />
        </Routes>
      </main>
      <Footer />
    </div>
  );
}
